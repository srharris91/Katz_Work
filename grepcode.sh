#!/bin/bash -i
clear
echo grep -rIn --exclude={*.old,*.svn,*.exec,*.vtu,*.pvtu,*.namelist*,*.shaun*,*SPTu*,*.C.*,*.dat,*.lay,*mesh} --exclude-dir={*.*,doc,Jon,tmp,Bumps} \"$1\" $*
grep -rIn --exclude={*.old,*.svn,*.exec,*.vtu,*.pvtu,*.namelist*,*.shaun*,*SPTu*,*.C.*,*.dat,*.lay,*.mesh} --exclude-dir={*.*,doc,Jon,tmp,Bumps} "$1" $*
