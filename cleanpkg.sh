#!/bin/bash
rm -f *~
rm -f codep/*~
rm -f codep/man/*~
rm -f codep/R/*~
rm -f codep/src/*~
rm -f codep/src/*.so
rm -f codep/src/*.o
rm -f codep/src/*.rds
cd codep && find -type f \( -not -name "MD5" \) -exec md5sum '{}' \; > MD5
