#!/bin/bash

CWPROOT=/opt/cwp/44R1/bin

rm -f *.ps *.eps *.png

$CWPROOT/psimage < snapz00801 \
    n1=481 d1=0.5 d1num=40 n1tic=2 f1=-20 f1min=-20 label1='Depth (m)' \
    n2=481 d2=0.5 d2num=40 n2tic=2 f2=-20 f2min=-20 label2='Distance (m)' \
    curve=pml.txt npair=5 curvedash=5 \
    clip=5.0e-11 d1s=0.5 d2s=0.5 \
    labelfont='SanSerif-Roman' labelsize=20 \
    height=4.0 width=4.0 > snapz00801_fs1.ps
ps2eps -l -B -s b0 -F snapz00801_fs1.ps
convert -density 300 snapz00801_fs1.eps snapz00801_fs1.png

$CWPROOT/psimage < snapz01201 \
    n1=481 d1=0.5 d1num=40 n1tic=2 f1=-20 f1min=-20 label1='Depth (m)' \
    n2=481 d2=0.5 d2num=40 n2tic=2 f2=-20 f2min=-20 label2='Distance (m)' \
    curve=pml.txt npair=5 curvedash=5 \
    clip=5.0e-11 d1s=0.5 d2s=0.5 \
    labelfont='SanSerif-Roman' labelsize=20 \
    height=4.0 width=4.0 > snapz01201_fs1.ps
ps2eps -l -B -s b0 -F snapz01201_fs1.ps
convert -density 300 snapz01201_fs1.eps snapz01201_fs1.png

$CWPROOT/psimage < snapz01601 \
    n1=481 d1=0.5 d1num=40 n1tic=2 f1=-20 f1min=-20 label1='Depth (m)' \
    n2=481 d2=0.5 d2num=40 n2tic=2 f2=-20 f2min=-20 label2='Distance (m)' \
    curve=pml.txt npair=5 curvedash=5 \
    clip=5.0e-11 d1s=0.5 d2s=0.5 \
    labelfont='SanSerif-Roman' labelsize=20 \
    height=4.0 width=4.0 > snapz01601_fs1.ps
ps2eps -l -B -s b0 -F snapz01601_fs1.ps
convert -density 300 snapz01601_fs1.eps snapz01601_fs1.png

rm -f *.ps *.eps
