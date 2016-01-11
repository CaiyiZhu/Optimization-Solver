function [n F]=find_n_FS(x)
fn_2=0;
fn_1=1;
fn=fn_1+fn_2;
n=2;
F=[fn_2;fn_1;fn];
while fn<x
    fn_2=fn_1;
    fn_1=fn;
    fn=fn_1+fn_2;
    n=n+1;
    F=[F;fn];
end