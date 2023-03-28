#!/bin/bash
function copyfile(){
	path=$(pwd)
	time=$(date "+%Y-%m-%d_%H:%M")
	Creat_path="${path}/${time}"
	echo ${Creat_path}
	mkdir ${Creat_path}
	cp $1 ${Creat_path}
	cp $2 ${Creat_path}
	cp $3 ${Creat_path}
	cp $4 ${Creat_path}
	cp $5 ${Creat_path}
	cd ${Creat_path}
	qsub run.sh
	shopt -s  expand_aliases 
	shopt expand_aliases
	alias draw='/home/mozhu/vampire-develop/Mn2Au/draw.sh'
	while true
	do
		if [ -e "log" ];then
			if (grep -q "Simulation ended" log);then
			curl "https://api2.pushdeer.com/message/push?pushkey=PDU20757TBbImxQwXVYSCQlRj7lD8JvOkUtVPf6t3&text=${Creat_path}_Simulation_finish"
			break
			elif (grep -q "Fatal error" log);then
			curl "https://api2.pushdeer.com/message/push?pushkey=PDU20757TBbImxQwXVYSCQlRj7lD8JvOkUtVPf6t3&text=${Creat_path}_Simulation_fail"
			break
			else
				continue
			fi
		else 
			continue
		fi
	done	

}

copyfile ./input ./vampire-serial ./Mn2Au.ucf ./Mn2Au.mat ./run.sh
