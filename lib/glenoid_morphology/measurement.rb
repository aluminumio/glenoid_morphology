# frozen_string_literal: true

module GlenoidMorphology
  # Result of `GlenoidMorphology.measure`.
  #
  # `pico_circle` is the inferior-2/3 reference circle (Pico method, Baudi 2005).
  # `bone_loss_percent` is computed from the active reference circle (Pico when
  # `method: :pico` -- the default -- otherwise the full-rim circle).
  # `legacy_bone_loss_percent` is the full-circle area-integral metric retained
  # for back-comparison; not part of the public stable API.
  Measurement = Struct.new(
    :en_face_plane,             # { a:, b:, c:, d: } -- ax + by + cz = d in mm space
    :glenoid_rim_points,        # Array<Array<Float, 3>> -- ordered 3D mm points on the rim
    :fitted_circle,             # { center:, radius_mm:, normal: } -- ACTIVE reference circle
    :pico_circle,               # { center:, radius_mm:, normal: } -- inferior-2/3 fit (nil if full_circle)
    :defect_arc,                # { start_angle:, end_angle:, arc_length_mm: }
    :bone_loss_percent,         # 0..100 (Pico-sector or full-circle, per method)
    :legacy_bone_loss_percent,  # 0..100 (full-circle area integral, for comparison)
    :method,                    # :pico or :full_circle
    :confidence,                # 0..1
    keyword_init: true
  )
end
