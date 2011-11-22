#!/bin/sh

root -b -q performFits_115.C | grep DATACARD | tr -d DATACARD > headSignal_M115.h
root -b -q performFits_110.C | grep DATACARD | tr -d DATACARD > headSignal_M110.h
root -b -q performFits_120.C | grep DATACARD | tr -d DATACARD > headSignal_M120.h
root -b -q performFits_130.C | grep DATACARD | tr -d DATACARD > headSignal_M130.h
root -b -q performFits_140.C | grep DATACARD | tr -d DATACARD > headSignal_M140.h
root -b -q performFits_Bg.C  | grep DATACARD | tr -d DATACARD > headBg.h
