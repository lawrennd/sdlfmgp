function addToolboxes(os,add)


if os==0,
    dir2 = '/kyb/agbs/malvarez/myfiles/mlprojects/';
else
   dir2 = 'D:\Algoritmos_Doctorado\mlprojects2\mlprojects\';
end

if add
    func = str2func('addpath');
else
    func = str2func('rmpath');
end

% Updated in SVN yet
func(strcat(dir2,'matlab/netlab/NETLAB3p3'))
func(strcat(dir2,'mocap/matlab'))
func(strcat(dir2,'mltools/matlab'))
func(strcat(dir2,'kern/matlab'))
func(strcat(dir2,'ndlutil/matlab'))
func(strcat(dir2,'gp/matlab'))
func(strcat(dir2,'optimi/matlab'))
func(strcat(dir2,'datasets/matlab'))
func(strcat(dir2,'multigp/matlab'))


