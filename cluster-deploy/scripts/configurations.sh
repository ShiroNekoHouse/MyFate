#!/bin/bash

user=app
dir=/data/projects/fate
mysqldir=/data/projects/common/mysql/mysql-8.0
javadir=/data/projects/common/jdk/jdk1.8
venvdir=/data/projects/fate/venv
redispass=fate1234
partylist=(partyA.id partyB.id) 
JDBC0=(A.MS-ip A.dbname A.user A.password) 
JDBC1=(B.MS-ip B.dbname B.user B.password) 
iplist=(A.F-ip A.MS-ip A.P-ip A.R-ip A.S-ip A.E-ip B.F-ip B.MS-ip B.P-ip B.R-ip B.S-ip B.E-ip)
fedlist0=(A.F-ip)
fedlist1=(B.F-ip)
meta0=(A.MS-ip)
meta1=(B.MS-ip)
proxy0=(A.P-ip)
proxy1=(B.P-ip)
roll0=(A.R-ip)
roll1=(B.R-ip)
egglist0=(A.E1-ip A.E2-ip A.E3-ip...)
egglist1=(B.E1-ip B.E2-ip B.E3-ip...) 
tmlist0=(A.TM-ip)
tmlist1=(B.TM-ip)
serving0=(A.S1-ip A.S2-ip)
serving1=(B.S1-ip B.S2-ip)
exchangeip=exchangeip 