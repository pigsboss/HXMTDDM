function hxmtassworker(inputfile,NUMIT)
%HXMTASSWORKER
if ischar(NUMIT)
  NUMIT=floor(str2double(NUMIT));
end

fileList=getFileList(inputfile);
NumFiles=length(fileList);
display([int2str(NumFiles),' jobs to process:'])
t=0;
for n=1:NumFiles
  tStart=tic;
  job=load(fileList{n});
  job.I=addmaz(job.PSF,job.obs,job.psi_obs,job.I,job.bg,[],NUMIT);
  job.t=job.t+NUMIT;
  t=toc(tStart)+t;
  display([int2str(n),' of ',int2str(NumFiles),...
    ' jobs done. Time elapsed: ',num2str(t),...
    ' seconds. Time remaining: ',num2str(t/n*(NumFiles-n)),' seconds...'])
  save(fileList{n},'-struct','job')
end
return

function fList=getFileList(inputfile)
%GETFILELIST
fid=fopen(inputfile,'r');

% 1st-pass test:
N=0;
tline=fgetl(fid);
while ischar(tline)
  if exist(tline,'file')==2
    N=N+1;
  end
  tline=fgetl(fid);
end

frewind(fid)
fList=cell(N,1);
tline=fgetl(fid);
n=1;
while ischar(tline)
  if exist(tline,'file')==2
    fList{n}=tline;
    n=n+1;
  end
  tline=fgetl(fid);
end
fclose(fid);
return