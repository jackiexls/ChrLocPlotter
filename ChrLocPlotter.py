#!/usr/bin/env python

'''
**--coding: utf-8--**
**--author: xinglongsheng@caas.cn--**
**--usage:  This script is used for drawing chromosomal location graph 
            for genes or other genomic features --**
**--date: 20190310--**
'''

from __future__ import division
import os
import sys
import re
import argparse
import svgwrite
import math


__version__ = 'V1.0'

helpInfo = '''
    This program is used to draw gene location plot on each chromosome based on the location of genes
of interest in the input file. The result for multiple chromosomes is integrated into single one canvas.
Meanwhile, this program also can be used for other types of genomic features, such as exon location, 
SNP position, molecular marks, etc. Importantly, the density of the features should not be too high so 
that they can be clearly displayed in the figure although some parameters could be adjusted to modify 
the appearance of the resultant figure.                                           
    This program receives two files as the input, namely the chromosome length file and gene chromosome 
location file. Both files should be provided in tab-delimited plain-text format.
'''

parser = argparse.ArgumentParser(description=helpInfo)
parser.add_argument('--length','-l',required=True,help='chromosome length file')
parser.add_argument('--position','-p',required=True,help='file containing gene position info')
parser.add_argument('--outfile','-o',required=True,help='prefix of output filename')
parser.add_argument('--rulerRatio','-r',required=True,help='a zooming ratio for drawing length scale')
argv = parser.parse_args()

mb = 10**6
kb = 10**3
cmRatio = 5


def readLocusFile(filename):
	ifh = open(filename)
	
	locus = {}
	while True:
		line = ifh.readline()
		if len(line) == 0:
			break
		line = line.rstrip()
		arr = line.split('\t')
		
		if not locus.has_key(arr[0]):
			locus[arr[0]] = []
		else:
			pass
		temp = {}

		temp['start'] = int(arr[1])
		temp['end'] = int(arr[2])
		temp['position'] = (temp['start'] + temp['end']) / 2.0 / mb * cmRatio if (temp['end'] - temp['start']) % 2 == 0 else (temp['start'] + temp['end'] - 1) / 2.0 / mb * cmRatio
		temp['name'] = arr[3]
		temp['strand'] = arr[4]
		locus[arr[0]].append(temp)

	ifh.close()
	return locus


def drawGenePosition(positionDict,outfile,chromLenDict,rulerRatio):
	xstart0 = -50
	xstart1 = -30
	ystart = 10
	topDist = 6.9
	bottomDist = 7.14
	rx = 5.15
	ry = 4.55
	columnNum = 7
	rowSpace = 5
	offset = 7
	rulerColor = "#808080"
	strokeColor = '#808080'
	fillColor = '#4F81BD'
	rulerWidth = '0.5'
	labelWidth = '0.5'
	chrWidth = 10.3
	fontsize = 5
	fontStyle = "font-size:%d;font-family:%s" % (fontsize, "Arial")

	dwg = svgwrite.Drawing(outfile,profile='full')
	xstart = xstart0
	rowNum = int(math.ceil(len(sorted(positionDict.keys())) / columnNum))
	tickRanges = []
	for i in range(0, rowNum):
		tickRanges.append(int(math.ceil(max([chromLenDict[x] for x in sorted(positionDict.iterkeys())[i*columnNum:(i+1)*columnNum]]) / cmRatio)))
	print(tickRanges)
	chromCount = 0
	for chrom in sorted(positionDict.iterkeys()):
		if chromCount % columnNum == 0:			
			print('Add the chromosome ruler on the left side')
			lastUpperRange = int(math.ceil(tickRanges[int(math.floor(chromCount / columnNum)) - 1] / 5)) * 5
			upperRange = int(math.ceil(tickRanges[int(math.floor(chromCount / columnNum))] / 5)) * 5
			print(chromCount, xstart, ystart, chromCount, lastUpperRange, upperRange)
			if chromCount > 0:
				ystart += topDist + (lastUpperRange + rowSpace) * cmRatio * rulerRatio
			dwg.add(dwg.line((20,ystart + topDist),(20,ystart + topDist + upperRange * cmRatio * rulerRatio),stroke=rulerColor,stroke_width=rulerWidth))
			g = dwg.g(style = fontStyle)
			g.add(dwg.text('Mb', insert = (15, ystart + topDist - 3)))
			dwg.add(g)

			for i in range(0, upperRange * cmRatio + 1, 1):
				if i % 5 == 0:
					dwg.add(dwg.line((15.7,ystart + topDist + i * rulerRatio),(20,ystart + topDist + i * rulerRatio),stroke=rulerColor,stroke_width=rulerWidth))
					g = dwg.g(style=fontStyle)
					g.add(dwg.text(round(i/cmRatio,1),insert=(2.7,ystart + topDist + i * rulerRatio + 2.5)))
					dwg.add(g)
				else:
					dwg.add(dwg.line((17.7,ystart + topDist + i * rulerRatio),(20,ystart + topDist + i * rulerRatio),stroke=rulerColor,stroke_width=rulerWidth))
			xstart = xstart1		
		tempPosArr = positionDict[chrom]
		positionArr = sorted(tempPosArr, key = lambda comp : comp['position'])
		xstart = xstart + 85
		height = chromLenDict[chrom]

		dwg.add(dwg.rect(insert = (xstart, ystart), size = (10.3,height * rulerRatio + topDist + bottomDist), rx = rx, ry = ry, stroke = '#000000', stroke_width='0.5',fill=fillColor,opacity='0.5'))
		g = dwg.g(style = fontStyle)
		g.add(dwg.text(chrom.capitalize(), insert = (xstart - 3, ystart - 3)))
		dwg.add(g)
		lx2 = xstart - 0.4
		lx1 = lx2 - 6.3	
		rx1 = xstart + chrWidth + 0.4
		rx2 = rx1 + 6.3
		
		dwg.add(dwg.line((lx1, ystart + topDist), (lx2, ystart + topDist), stroke = strokeColor, stroke_width = labelWidth))
		
		g = dwg.g(style=fontStyle)
		g.add(dwg.text('0',insert=(xstart - 14,ystart + topDist + 3)))
		dwg.add(g)
		flagList = []
		flagNum = 0
		lastYpos = ystart + topDist + 3

		for i in range(0,len(positionArr)):
			y1 = y2 = ystart + topDist + positionArr[i]['position'] * rulerRatio
			if i > 0 and positionArr[i]['position'] - positionArr[i-1]['position'] < 1.0:
				if y1 > lastYpos:
					dwg.add(dwg.line((lx1,y1),(lx2,y2),stroke=strokeColor,stroke_width=labelWidth))
					g = dwg.g(style=fontStyle)
					g.add(dwg.text(round(positionArr[i]['position']/cmRatio,2),insert=(lx1 - 22,y1 + 3)))
					dwg.add(g)

					dwg.add(dwg.line((rx1,y1),(rx2,y2),stroke=strokeColor,stroke_width=labelWidth))
					g = dwg.g(style=fontStyle)
					g.add(dwg.text(positionArr[i]['name'],insert=(rx2 + 2,y1 + 3)))
					dwg.add(g)
					lastYpos = y1 + 3
				else:
					dwg.add(dwg.polyline(points=((lx1,lastYpos + offset),(lx1 + 2.8,lastYpos + offset),(lx2 - 2.3,y1),(lx2,y1)),fill='none',stroke=strokeColor,stroke_width=labelWidth))
					g = dwg.g(style=fontStyle)
					g.add(dwg.text(round(positionArr[i]['position']/cmRatio,2),insert=(lx1 - 22, lastYpos + offset)))
					dwg.add(g)
					dwg.add(dwg.polyline(points=((rx1,y2),(rx1 + 2.3, y2),(rx2 - 2.8,lastYpos + offset),(rx2,lastYpos + offset)),fill='none',stroke=strokeColor,stroke_width=labelWidth))
					g = dwg.g(style=fontStyle)
					g.add(dwg.text(positionArr[i]['name'],insert=(rx2 + 2,lastYpos + offset)))
					dwg.add(g)
					lastYpos = lastYpos + offset
			else:
				flagNum = 0
				if y1 > lastYpos:
					dwg.add(dwg.line((lx1,y1),(lx2,y2),stroke=strokeColor,stroke_width=labelWidth))						
					g = dwg.g(style=fontStyle)
					g.add(dwg.text(round(positionArr[i]['position']/cmRatio,2),insert=(lx1 - 22,y1 + 3)))
					dwg.add(g)
					dwg.add(dwg.line((rx1,y1),(rx2,y2),stroke=strokeColor,stroke_width=labelWidth))
					g = dwg.g(style=fontStyle)
					g.add(dwg.text(positionArr[i]['name'],insert=(rx2 + 2,y1 + 3)))
					dwg.add(g)
					lastYpos = y1 + 3
				else:
					dwg.add(dwg.polyline(points=((lx1,lastYpos + offset),(lx1 + 2.8,lastYpos + offset),(lx2 - 2.3,y1),(lx2,y1)),fill='none',stroke=strokeColor,stroke_width=labelWidth))
					g = dwg.g(style=fontStyle)
					g.add(dwg.text(round(positionArr[i]['position']/cmRatio,2),insert=(lx1 - 22, lastYpos + offset)))
					dwg.add(g)
					dwg.add(dwg.polyline(points=((rx1,y2),(rx1 + 2.3,y2),(rx2 - 2.8,lastYpos + offset),(rx2,lastYpos + offset)),fill='none',stroke=strokeColor,stroke_width=labelWidth))
					g = dwg.g(style=fontStyle)
					g.add(dwg.text(positionArr[i]['name'],insert=(rx2 + 2, lastYpos + offset)))
					dwg.add(g)
					lastYpos = lastYpos + offset
		
		yend = ystart + height * rulerRatio + topDist
		if height - positionArr[i]['position'] < 1.0:
			dwg.add(dwg.polyline(points=((lx1,yend + offset),(lx1 + 2.8,yend + offset),(lx2 - 2.3,yend),(lx2,yend)),fill='none',stroke=strokeColor,stroke_width=labelWidth))
			g = dwg.g(style = fontStyle)
			g.add(dwg.text(round(height/cmRatio,2),insert=(lx1 - 22,yend + 10)))
			dwg.add(g)
		else:
			dwg.add(dwg.line((lx1,yend),(lx2,yend),stroke=strokeColor,stroke_width=labelWidth))
			if len(str(round(height,2))) <= 5:
				g = dwg.g(style = fontStyle)
				g.add(dwg.text(round(height/cmRatio,2),insert=(lx1 - 22,yend + 3)))
				dwg.add(g)
			else:
				g = dwg.g(style = fontStyle)
				g.add(dwg.text(round(height/cmRatio,2),insert=(lx1 - 22,yend + 3)))
				dwg.add(g)
		chromCount += 1
		
	dwg.save()


def readChromLength(filename):
	ifh = open(filename)
	chromDict = {}

	while True:
		line = ifh.readline()
		line = line.rstrip()
		if len(line) == 0:
			break
		arr = line.split('\t')
		chromDict[arr[0]] = int(arr[1]) / mb * cmRatio
	ifh.close()

	return chromDict


if __name__ == '__main__':
	chrlen = argv.length
	posfile = argv.position
	outfile = argv.outfile
	rulerRatio = float(argv.rulerRatio)

	print('Read in chromosome length info')
	chromLenDict = readChromLength(chrlen)
	print('Read in SNP position info')
	posDict = readLocusFile(posfile)
	outputfile = outfile + '.svg'
	print('Drawing gene position start')
	drawGenePosition(posDict, outputfile, chromLenDict, rulerRatio)
	print('Result output to %s' % outputfile)
	print('Drawing gene position end')
