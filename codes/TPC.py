#!/usr/bin/env python
#_*_coding:utf-8_*_

import re

def TPC(fastas, **kw):
	triPeptides = kw['AAA']
	encodings = []
	header = ['#'] + triPeptides
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		tmpCode = [0] * len(triPeptides)
		for j in range(len(triPeptides)):
			tmpCode[j] = sequence.count(triPeptides[j])
		if sum(tmpCode) != 0:
			tmpCode = [i/(len(sequence) - 2) for i in tmpCode]
		code = code + tmpCode
		encodings.append(code)
	return encodings
