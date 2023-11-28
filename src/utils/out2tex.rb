#!/usr/bin/env ruby
PREFIX="output"

File.open("#{PREFIX}/G.out") do |f|
  str = f.read.split("\n").map{|s| s.split(/\s+/)};
  File.open("#{PREFIX}/G.tex", "w") do |of|
    str.each_with_index{|row, id| of.write(row.join(" & ") + "  \\\\ \n")}
  end
end

