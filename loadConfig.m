function loadConfig
%LOADCONFIG Load the config file content of HXMT DDM.

global config;

if isempty(config)
  if exist('./config.mat','file')==2
    config=load('./config.mat');
  else
    disp('config.mat does not exist. Run configHXMT.m to build it.')
    config=configHXMTDDM;
  end
end
return
