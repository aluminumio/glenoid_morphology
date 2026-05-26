# frozen_string_literal: true

module GlenoidMorphology
  # Result of `GlenoidMorphology.measure`.
  Measurement = Struct.new(
    :en_face_plane,       # { a:, b:, c:, d: } -- ax + by + cz = d in mm space
    :glenoid_rim_points,  # Array<Array<Float, 3>> -- ordered 3D mm points on the rim
    :fitted_circle,       # { center:, radius_mm:, normal: }
    :defect_arc,          # { start_angle:, end_angle:, arc_length_mm: }
    :bone_loss_percent,   # 0..100
    :confidence,          # 0..1
    keyword_init: true
  )
end
