function [real_matrix] = complex2real(complex_matrix)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    real_matrix = [real(complex_matrix),-imag(complex_matrix);imag(complex_matrix),real(complex_matrix)];
end