@echo off
REM This part of batch file added by SVDAPREP
REM
REM Delete model input files.
del p2p_bigelk.dat > nul
REM
REM Run PARCALC to compute base parameters from super parameters.
parcalc > nul
REM
REM The following is copied directly from file model.bat
REM
@echo off
rem del %~dp0model.out
par2par p2p_bigelk.dat > nul
call m:\models\bacteria\hspf\bigelkhydrocal201601\hspf-files\winhspf.bat > nul
call m:\models\bacteria\hspf\bigelkhydrocal201601\r-files\post-proc.bat > nul
