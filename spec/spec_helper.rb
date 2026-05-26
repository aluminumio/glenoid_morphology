# frozen_string_literal: true

require "glenoid_morphology"
require_relative "support/synthetic"

RSpec.configure do |config|
  config.expect_with(:rspec) { |c| c.syntax = :expect }
  config.disable_monkey_patching!
  config.example_status_persistence_file_path = ".rspec_status"
  config.include Synthetic
end
