# frozen_string_literal: true

require "numo/narray"

# Synthetic scapula/glenoid generator for specs.
#
# A "scapula" is approximated as:
#   - A thin rectangular slab (the scapula body) along the YZ plane
#   - A disc-shaped cap at one end of the slab (the glenoid). The disc has
#     a chosen radius, lies in a plane with chosen normal, and has thickness
#     `glenoid_thickness` voxels.
#
# Test variants:
#   - intact_disc
#   - chord_defect(angle_deg:) — removes the half-plane on one side of a chord
#     whose subtended angle at the disc center equals `angle_deg`. This is the
#     classic geometric model of anterior glenoid bone loss (Hovelius / Friedman).
#
# All masks are returned as Numo::Bit shape [D, H, W].
module Synthetic
  module_function

  GRID = 100
  GLENOID_RADIUS = 24.0      # mm == voxels (we use 1mm spacing)
  GLENOID_THICKNESS = 4      # voxels
  BODY_THICKNESS = 2         # voxels (super-thin "scapular blade")
  BODY_HEIGHT = 60           # voxels
  GLENOID_CENTER = [50, 50, 80].freeze  # i, j, k (k is the lateral axis)

  # @return [Numo::Bit]
  def intact_glenoid_mask(grid: GRID,
                          radius: GLENOID_RADIUS,
                          thickness: GLENOID_THICKNESS,
                          center: GLENOID_CENTER,
                          plane_normal: [0.0, 0.0, 1.0])
    build_mask(grid: grid, radius: radius, thickness: thickness,
               center: center, plane_normal: plane_normal,
               defect_angle_deg: 0.0)
  end

  # Defect angle measured at the disc centre (subtended angle of the chord).
  def chord_defect_mask(angle_deg:,
                        grid: GRID,
                        radius: GLENOID_RADIUS,
                        thickness: GLENOID_THICKNESS,
                        center: GLENOID_CENTER,
                        plane_normal: [0.0, 0.0, 1.0],
                        defect_direction: [-1.0, 0.0, 0.0])
    build_mask(grid: grid, radius: radius, thickness: thickness,
               center: center, plane_normal: plane_normal,
               defect_angle_deg: angle_deg,
               defect_direction: defect_direction)
  end

  # Geometric expected bone-loss percent for a chord subtending `angle_deg`.
  # The chord is at distance d = r*cos(angle/2) from the centre; the segment area is
  #   A = 0.5 * r^2 * (theta - sin theta), theta = angle in rad.
  # As fraction of pi*r^2 it depends ONLY on theta:
  #   loss = (theta - sin theta) / (2*pi).
  def expected_loss_fraction(angle_deg)
    theta = angle_deg * Math::PI / 180.0
    (theta - Math.sin(theta)) / (2.0 * Math::PI)
  end

  # ---- internals ----

  def build_mask(grid:, radius:, thickness:, center:, plane_normal:,
                 defect_angle_deg:, defect_direction: [-1.0, 0.0, 0.0])
    mask = Numo::Bit.zeros(grid, grid, grid)

    n = normalise(plane_normal)
    d_dir = project_to_plane(defect_direction, n)
    d_dir = normalise(d_dir)
    chord_offset = (defect_angle_deg.positive? ? radius * Math.cos(defect_angle_deg * Math::PI / 360.0) : nil)

    grid.times do |i|
      grid.times do |j|
        grid.times do |k|
          dx = i - center[0]
          dy = j - center[1]
          dz = k - center[2]

          # Distance along plane normal -> thickness check
          along_n = dx * n[0] + dy * n[1] + dz * n[2]
          next if along_n.abs > (thickness / 2.0)

          # Distance from centre, in the plane
          in_plane_x = dx - along_n * n[0]
          in_plane_y = dy - along_n * n[1]
          in_plane_z = dz - along_n * n[2]
          rad = Math.sqrt(in_plane_x**2 + in_plane_y**2 + in_plane_z**2)
          next if rad > radius

          # Chord defect: cut off the half-plane d_dir > chord_offset.
          if chord_offset
            along_defect = in_plane_x * d_dir[0] + in_plane_y * d_dir[1] + in_plane_z * d_dir[2]
            next if along_defect > chord_offset
          end

          mask[i, j, k] = 1
        end
      end
    end

    mask
  end

  # Pear-shaped glenoid: a CIRCLE of radius `r` for the inferior 2/3 of the
  # rim, with the superior 1/3 narrowed inward by `taper_mm` to model the
  # natural pear shape of a healthy human glenoid.
  #
  # The in-plane basis (u, v) is picked by `in_plane_basis(plane_normal)`.
  # The SI axis is `v`: positive v = superior, negative v = inferior. The
  # superior narrowing kicks in for points with v > r/3 (i.e. the superior
  # third), and tapers the rim from radius `r` (at v = r/3) down to
  # `r - taper_mm` (at v = r, the superior pole).
  #
  # If you fit a circle to ALL the rim points of this shape you'll get a
  # radius somewhere between (r - taper) and r — that's the legacy failure
  # mode. The Pico method, by restricting to the inferior 240°, recovers
  # the true r.
  def pear_glenoid_mask(grid: GRID,
                        r: 12.0,
                        taper_mm: 4.0,
                        thickness: GLENOID_THICKNESS,
                        center: GLENOID_CENTER,
                        plane_normal: [0.0, 0.0, 1.0],
                        defect_angle_deg: 0.0,
                        defect_at_inferior: true)
    mask = Numo::Bit.zeros(grid, grid, grid)
    n = normalise(plane_normal)
    u_ax, v_ax = in_plane_basis(n)
    taper_start_v = r / 3.0
    taper_span    = r - taper_start_v # the v-range over which we narrow

    grid.times do |i|
      grid.times do |j|
        grid.times do |k|
          dx = i - center[0]; dy = j - center[1]; dz = k - center[2]
          along_n = dx * n[0] + dy * n[1] + dz * n[2]
          next if along_n.abs > (thickness / 2.0)

          in_plane_x = dx - along_n * n[0]
          in_plane_y = dy - along_n * n[1]
          in_plane_z = dz - along_n * n[2]
          u_pos = in_plane_x * u_ax[0] + in_plane_y * u_ax[1] + in_plane_z * u_ax[2]
          v_pos = in_plane_x * v_ax[0] + in_plane_y * v_ax[1] + in_plane_z * v_ax[2]

          # Below the superior taper start, the shape is the circle of radius r.
          rad2 = u_pos * u_pos + v_pos * v_pos
          next if rad2 > r * r

          # Above the taper start, narrow the |u| extent linearly toward the
          # superior pole.
          if v_pos > taper_start_v
            t = ((v_pos - taper_start_v) / taper_span).clamp(0.0, 1.0)
            u_limit = Math.sqrt([r * r - v_pos * v_pos, 0.0].max) - taper_mm * t
            next if u_pos.abs > [u_limit, 0.0].max
          end

          # Optional chord defect at the inferior pole (or superior pole).
          # The chord is at distance r * cos(angle/2) from the centre.
          if defect_angle_deg.positive?
            cut_distance = r * Math.cos(defect_angle_deg * Math::PI / 360.0)
            along = defect_at_inferior ? -v_pos : v_pos
            next if along > cut_distance
          end

          mask[i, j, k] = 1
        end
      end
    end
    mask
  end

  # Synthetic humerus mask: a sphere sitting in the +plane_normal direction
  # past the glenoid disc. Used to disambiguate which side of the disc faces
  # the joint. Returned as a Numo::Bit of the same shape.
  def humerus_mask_for(grid: GRID,
                       glenoid_center: GLENOID_CENTER,
                       plane_normal: [0.0, 0.0, 1.0],
                       offset: 12.0,
                       radius: 14.0)
    n = normalise(plane_normal)
    c = [
      glenoid_center[0] + offset * n[0],
      glenoid_center[1] + offset * n[1],
      glenoid_center[2] + offset * n[2]
    ]
    mask = Numo::Bit.zeros(grid, grid, grid)
    r2 = radius * radius
    grid.times do |i|
      grid.times do |j|
        grid.times do |k|
          dx = i - c[0]; dy = j - c[1]; dz = k - c[2]
          mask[i, j, k] = 1 if dx * dx + dy * dy + dz * dz <= r2
        end
      end
    end
    mask
  end

  def normalise(v)
    n = Math.sqrt(v.inject(0.0) { |s, x| s + x * x })
    v.map { |x| x / n }
  end

  def project_to_plane(v, n)
    dot = v[0] * n[0] + v[1] * n[1] + v[2] * n[2]
    [v[0] - dot * n[0], v[1] - dot * n[1], v[2] - dot * n[2]]
  end

  def in_plane_basis(n)
    seed = (n[0].abs > 0.9 ? [0.0, 1.0, 0.0] : [1.0, 0.0, 0.0])
    u = project_to_plane(seed, n)
    u = normalise(u)
    v = cross(n, u)
    [u, normalise(v)]
  end

  def cross(a, b)
    [a[1] * b[2] - a[2] * b[1],
     a[2] * b[0] - a[0] * b[2],
     a[0] * b[1] - a[1] * b[0]]
  end
end
