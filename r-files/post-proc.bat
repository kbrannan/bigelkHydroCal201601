echo off
del model.ins
del model.out
"C:\Program Files\R\R-3.1.3\bin\x64\Rscript.exe" --vanilla %~dp0hspf-output-proc.R