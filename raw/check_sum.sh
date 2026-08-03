#!/bin/bash
xargs -P 36 -d '\n' -I{} bash -c 'echo "{}" | md5sum -c -' < md5.txt
