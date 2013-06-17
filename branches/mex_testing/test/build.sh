#!/bin/sh
cd ../ && make -j2 install && cd test && make sttsm_but_one
