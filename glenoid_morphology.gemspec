# frozen_string_literal: true

require_relative "lib/glenoid_morphology/version"

Gem::Specification.new do |spec|
  spec.name        = "glenoid_morphology"
  spec.version     = GlenoidMorphology::VERSION
  spec.authors     = ["Jonathan Siegel"]
  spec.email       = ["jonathan@siegel.io"]

  spec.summary     = "Geometric glenoid bone-loss measurement from shoulder CT segmentations"
  spec.description = "Pure-Ruby implementation of the en-face / best-fit-circle method for " \
                     "measuring anterior glenoid bone loss from a scapula label volume. " \
                     "Ports the geometry pipeline described in arxiv 2511.14083. " \
                     "Consumes the per-voxel masks produced by `shoulder_segmenter`."
  spec.homepage    = "https://github.com/aluminumio/glenoid_morphology"
  spec.license     = "MIT"

  spec.required_ruby_version = ">= 3.3.0"

  spec.metadata["homepage_uri"]    = spec.homepage
  spec.metadata["source_code_uri"] = spec.homepage
  spec.metadata["bug_tracker_uri"] = "#{spec.homepage}/issues"

  spec.files = Dir[
    "lib/**/*.rb",
    "doc/*.md",
    "README.md",
    "LICENSE",
    "glenoid_morphology.gemspec"
  ]
  spec.require_paths = ["lib"]

  spec.add_dependency "numo-narray", "~> 0.9"
  spec.add_dependency "numo-linalg", "~> 0.1"

  spec.add_development_dependency "bundler", ">= 2.0"
  spec.add_development_dependency "rspec",   "~> 3.13"
end
