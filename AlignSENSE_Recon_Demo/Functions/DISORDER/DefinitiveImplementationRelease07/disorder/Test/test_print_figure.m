%========================================================================================= 
% Program: print_figure.m % 
% Usage: % matlab -nodesktop -nodisplay -nosplash -r "test_print_figure('file_name','file_format');exit" 
%========================================================================================= 
function [] = test_print_figure( outfile, file_format ) 
disp('A simple test to illustrate generating figures in a batch mode.'); 
x = 0:.1:1; 
A = exp(x); 
plot(x,A); 
print(outfile,file_format); 
end