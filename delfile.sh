#!/bin/bash
delfile(){
path=$(pwd)
for i in `ls -d */`
do
	cd ${path}/$i
	if (grep -q "Fatal error" log);then
		cd ..
		echo $i
		rm -r $i
	elif [ ! -e "log" ];then
		cd ..
		echo $i
		rm -r $i
	fi
done

}

delfile
