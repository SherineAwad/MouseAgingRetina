#!/bin/bash

 awk '$6 == "exonic" {print $7}' 106w_gex_hc_somatic.mm10_multianno.txt > W106_gex_exonic_hc_somatic.txt 
 grep -Ff W106_gex_exonic_hc_somatic.txt allDEGs.csv >> W106_gex_exonic_somatic_nearby.txt

 awk '$6 == "exonic" {print $7}' 108w_gex_hc_somatic.mm10_multianno.txt > W108_gex_exonic_hc_somatic.txt
 grep -Ff W108_gex_exonic_hc_somatic.txt allDEGs.csv >> W108_gex_exonic_somatic_nearby.txt
