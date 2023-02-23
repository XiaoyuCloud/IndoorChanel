function [ ODELTA ] = delta( t,t0 )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
if t>15&&t<16
    if t==15.4
        ODELTA=1;
    else
        ODELTA=0;
    end;
elseif abs(t-t0)<=0.1
    ODELTA=1;
else
    ODELTA=0;
end;
    
end

