#!/bin/bash
echo "NanoShaper from Docker using VMD"

export CONF=$1
AAA=$(grep Root_FileName $1)
echo $AAA
set -- $AAA
DIR=$(dirname "$3")
echo $DIR
docker run -it --mount type=bind,source="$DIR",target="$DIR" nanoshaper $CONF