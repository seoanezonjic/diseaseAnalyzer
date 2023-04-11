#! /usr/bin/env ruby

require 'optparse'
require 'json'

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:input] = nil
  opts.on("-i", "--input_file PATH", "Path to input file") do |item|
    options[:input] = item
  end

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end
end.parse!

data = JSON.parse(File.open(options[:input]).read)
items = []
data.dig('associations').each do |assoc|
	 items << assoc['object']['id']
end


items.each do |item|
	puts item
end