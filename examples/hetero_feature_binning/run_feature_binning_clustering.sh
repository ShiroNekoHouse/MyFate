#!/usr/bin/env bash
#
#  Copyright 2019 The FATE Authors. All Rights Reserved.
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#

cur_dir=$(pwd)

role=${1}
guest_id=${2}
host_id=${3}
jobid=${4}

nohup python ${cur_dir}/run_binning.py 1 ${jobid} ${role} ${guest_id} ${host_id} > nohup.${role} 2>&1 &

log_path=${cur_dir}/../../logs/${jobid}
echo "Please check log in " ${log_path}