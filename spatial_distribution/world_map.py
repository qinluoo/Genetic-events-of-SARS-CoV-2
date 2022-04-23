# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 14:12:51 2021

@author: qinluyao
"""


from pyecharts import Map
from sys import argv
f = open(argv[1], "r")
lines = f.readlines()
value = []
attr = []
for country in lines[0].strip().split("\t"):
    attr.append(country)
for count in lines[1].strip().split("\t"):
    value.append(int(count))
map0 = Map("", width=1200, height=600)
map0.add("Recombinant of " + argv[1][4:-4], attr, value, visual_range=[0, max(value)], maptype="world", is_visualmap=True, visual_text_color='#000', is_map_symbol_show=False)
map0.render(path=".".join(argv[1].split(".")[:-1]) + ".html")