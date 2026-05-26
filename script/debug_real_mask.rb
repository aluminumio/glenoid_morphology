#!/usr/bin/env ruby
# frozen_string_literal: true
#
# Instruments GlenoidMorphology.measure on a real shoulder segmentation,
# stage-by-stage. Used to diagnose the v0.1.0 -> v0.2.0 fix (whole-scapula
# rim fit on real masks).
#
# Usage:
#   source ~/.rvm/environments/ext-ruby-3.3.11
#   bundle exec ruby script/debug_real_mask.rb
#
# Inputs (hard-coded so the failure is reproducible):
#   /Users/jonathan/Projects/dicom-imager/tmp/segmentation_run_appendicular/labels.nii.gz
#   (scapula_left = label id 3, humerus_left = 1, femur_left = 7 in shoulder_segmenter task 294)

$LOAD_PATH.unshift(File.expand_path("../lib", __dir__))
$LOAD_PATH.unshift("/Users/jonathan/Projects/nifti-ruby/lib")

require "json"
require "fileutils"
require "numo/narray"
require "nifti"
require "glenoid_morphology"

LABELS_PATH = "/Users/jonathan/Projects/dicom-imager/tmp/segmentation_run_appendicular/labels.nii.gz"
OUT_DIR     = File.expand_path("debug_out", __dir__)
FileUtils.mkdir_p(OUT_DIR)

def log(msg) = warn(msg)

def write_csv(path, rows)
  File.open(path, "w") do |f|
    rows.each { |r| f.puts(r.map(&:to_s).join(",")) }
  end
end

def bbox(idx)
  return nil if idx.shape[0].zero?

  v = Numo::DFloat.cast(idx)
  mn = v.min(axis: 0).to_a
  mx = v.max(axis: 0).to_a
  { min: mn, max: mx, span: mx.zip(mn).map { |a, b| a - b } }
end

# --- Load NIfTI labels ---
log "Loading labels: #{LABELS_PATH}"
vol = Nifti.load(LABELS_PATH)
log "  shape=#{vol.shape} dtype=#{vol.dtype} voxel_size=#{vol.voxel_size.inspect}"
log "  affine=#{vol.affine.inspect}"

flat = vol.to_a
sx, sy, sz = vol.shape
nx, ny, nz = sx, sy, sz
log "  unpacking #{flat.size} voxels into (#{nx},#{ny},#{nz}) (x,y,z) on-disk order"

data_zyx = Numo::UInt8.cast(flat).reshape(nz, ny, nx)
log "  data_zyx shape=#{data_zyx.shape} (this is z,y,x)"

# The pipeline (run_segmentation.rb) passes seg.mask_for, where seg.data is
# in (z,y,x) order and seg.affine is in (x,y,z). The mismatch matters for
# physical-mm calculations but for visual stage analysis we'll use (z,y,x)
# throughout — that's exactly what the failing pipeline saw.

LABELS = {
  "scapula_left"  => 3,
  "humerus_left"  => 1,
  "femur_left"    => 7
}.freeze

masks = {}
LABELS.each do |name, id|
  m = data_zyx.eq(id).cast_to(Numo::Bit)
  cnt = m.count_true
  log "label=#{name} id=#{id} voxels=#{cnt}"
  masks[name] = m
end

hum_union = (masks["humerus_left"] | masks["femur_left"]).cast_to(Numo::Bit)
log "humerus_union (humerus_left|femur_left) voxels=#{hum_union.count_true}"

scap = masks["scapula_left"]
affine = vol.affine

# --- Largest connected component via iterative DFS ---
def largest_cc(mask)
  bit = mask.cast_to(Numo::Bit)
  shape = bit.shape
  flat_true = bit.where.to_a
  in_set = {}
  flat_true.each { |fi| in_set[fi] = true }
  visited = {}
  hw = shape[1] * shape[2]
  w  = shape[2]
  best_size = 0
  best_set = nil
  flat_true.each do |start|
    next if visited[start]

    stack = [start]
    visited[start] = true
    component = []
    while (cur = stack.pop)
      component << cur
      i = cur / hw
      r = cur % hw
      j = r / w
      k = r % w
      [[i + 1, j, k], [i - 1, j, k], [i, j + 1, k], [i, j - 1, k], [i, j, k + 1], [i, j, k - 1]].each do |ni, nj, nk|
        next if ni < 0 || nj < 0 || nk < 0
        next if ni >= shape[0] || nj >= shape[1] || nk >= shape[2]

        fi = ni * hw + nj * w + nk
        next if visited[fi]
        next unless in_set[fi]

        visited[fi] = true
        stack << fi
      end
    end
    if component.size > best_size
      best_size = component.size
      best_set = component
    end
  end
  [best_size, best_set, shape]
end

log ""
log "=== stage 1b: connected components ==="
t0 = Time.now
best_size, best_set, shape = largest_cc(scap)
log "scapula_left total=#{scap.count_true}, largest CC=#{best_size} (#{(best_size.to_f / scap.count_true * 100).round(2)}%) in #{(Time.now - t0).round(1)}s"

cc_mask = Numo::Bit.zeros(*shape)
hw = shape[1] * shape[2]; w = shape[2]
best_set.each do |fi|
  i = fi / hw; rem = fi % hw; j = rem / w; k = rem % w
  cc_mask[i, j, k] = 1
end

# --- Spatial relationship between scapula CC and humerus ---
log ""
log "=== stage 2b: scapula vs humerus geometry ==="
scap_idx = GlenoidMorphology::SurfaceExtraction.indices_of_true(cc_mask)
hum_idx  = GlenoidMorphology::SurfaceExtraction.indices_of_true(hum_union)
scap_c   = Numo::DFloat.cast(scap_idx).mean(axis: 0)
hum_c    = Numo::DFloat.cast(hum_idx).mean(axis: 0)
log "  scapula CC centroid (voxel)= #{scap_c.to_a.map { |x| x.round(2) }.inspect}"
log "  humerus_union centroid    = #{hum_c.to_a.map { |x| x.round(2) }.inspect}"
log "  vec (hum - scap, voxels)  = #{(hum_c - scap_c).to_a.map { |x| x.round(2) }.inspect}"
log "  scapula bbox: #{bbox(scap_idx).inspect}"
log "  humerus bbox: #{bbox(hum_idx).inspect}"

# --- Find voxels CLOSE TO humerus: candidate glenoid region ---
log ""
log "=== stage 2c: voxels close to humerus centroid (proximity heuristic) ==="

# For every surface voxel of the scapula, compute distance to humerus centroid.
# Then take the closest 5%, 10%, 20% and see how compact they are.
scap_surface = GlenoidMorphology::SurfaceExtraction.surface_voxels(cc_mask)
log "  scapula CC surface voxels=#{scap_surface.shape[0]}"

# Distance from each scap-surface voxel to the humerus CENTROID, in mm
# (using the affine).
affine_n = GlenoidMorphology::Affine.to_narray(affine)
voxel_spacing = GlenoidMorphology::Affine.voxel_spacing_mm(affine_n)
log "  estimated isotropic spacing ≈ #{voxel_spacing.round(3)} mm/voxel"

scap_surface_mm = GlenoidMorphology::Affine.voxels_to_mm(scap_surface, affine)
# But hum_c is voxel; bring it to mm too.
hum_c_mm_full = affine_n.dot(Numo::DFloat[hum_c[0], hum_c[1], hum_c[2], 1.0])
hum_c_mm = Numo::DFloat[hum_c_mm_full[0], hum_c_mm_full[1], hum_c_mm_full[2]]
log "  humerus centroid (mm) = #{hum_c_mm.to_a.map { |x| x.round(2) }.inspect}"

# Distances
dx = scap_surface_mm[true, 0] - hum_c_mm[0]
dy = scap_surface_mm[true, 1] - hum_c_mm[1]
dz = scap_surface_mm[true, 2] - hum_c_mm[2]
dist = Numo::NMath.sqrt(dx * dx + dy * dy + dz * dz)
log "  dist to humerus centroid: min=#{dist.min.round(2)} max=#{dist.max.round(2)} mean=#{dist.mean.round(2)} mm"

sorted = dist.sort
[0.02, 0.05, 0.10, 0.20].each do |q|
  cutoff = sorted[(sorted.size * q).to_i]
  mask = dist <= cutoff
  idx = mask.where
  kept = scap_surface[idx, true]
  log "  closest #{(q * 100).to_i}% (n=#{kept.shape[0]}, dist<=#{cutoff.round(2)}mm) bbox(vox)=#{bbox(kept).inspect}"
  kept_mm = scap_surface_mm[idx, true]
  cmm = kept_mm.mean(axis: 0).to_a.map { |x| x.round(2) }
  smm = kept_mm.max(axis: 0).to_a.zip(kept_mm.min(axis: 0).to_a).map { |a, b| (a - b).round(2) }
  log "    mm centroid=#{cmm.inspect}, mm span=#{smm.inspect}"
end

# --- Run full measure on each variant and report ---
log ""
log "=== stage 9: end-to-end measure ==="

def run_one(label, mask, humerus:, affine:, side:)
  t0 = Time.now
  result = GlenoidMorphology.measure(
    scapula_mask: mask, humerus_mask: humerus,
    affine: affine, side: side, seed: 12345
  )
  elapsed = (Time.now - t0).round(2)
  warn "  [#{label}]  blp=#{result.bone_loss_percent.round(2)}%  r=#{result.fitted_circle[:radius_mm].round(2)}mm  conf=#{result.confidence.round(3)}  (#{elapsed}s)"
  result
rescue => e
  warn "  [#{label}]  FAILED: #{e.class}: #{e.message}"
  warn "    #{e.backtrace.first(3).join("\n    ")}"
  nil
end

run_one "raw_scapula  / no humerus",    scap,    humerus: nil,       affine: affine, side: :left
run_one "raw_scapula  / hum_union",     scap,    humerus: hum_union, affine: affine, side: :left
run_one "largestCC    / no humerus",    cc_mask, humerus: nil,       affine: affine, side: :left
run_one "largestCC    / hum_union",     cc_mask, humerus: hum_union, affine: affine, side: :left

log ""
log "done"
