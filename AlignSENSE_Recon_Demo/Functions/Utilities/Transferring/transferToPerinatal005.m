

function [] = transferToPerinatal005(dirIn, dirOu)

command = sprintf('scp -r /home/ybr19/%s ybr19@perinatal005-pc:/home/ybr19/%s',dirIn,dirOu);
dos(command);

