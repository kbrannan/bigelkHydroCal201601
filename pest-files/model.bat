@echo off
rem del %~dp0model.out
par2par p2p_bigelk.dat > nul
call m:\models\bacteria\hspf\bigelkhydrocal201601\hspf-files\winhspf.bat > nul
call m:\models\bacteria\hspf\bigelkhydrocal201601\r-files\post-proc.bat > nul