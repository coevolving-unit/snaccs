#!/bin/bash

wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install
perl configureHomer.pl -install hg38
cp  HOCOMOCOv11_full_HUMAN_mono_homer_format_0.001.motif motifs/HOCOMOCOv11.motif
cp  HOCOMOCOv11_full_HUMAN_mono_homer_format_0.001.motif data/knownTFs/HOCOMOCOv11/all.motifs
