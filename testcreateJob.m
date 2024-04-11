clc
clear
% s = toolboxdir('simulink');
% rmpath(s);
profile on
myCluster = parcluster("Processes");
j = createJob(myCluster); 
N = 2000;
createTask(j, @A, 1, {N});
createTask(j, @B, 1, {2*N});
submit(j);
wait(j);
result_1 = fetchOutputs(j);
destroy(j);

A(N);
B(2*N);

profile viewer

function state = A(input)
    state = 0;
    for i = 1:input
        for j = 1:input
        state = state + i*1e-6;
        end
    end
end

function state = B(input)
    state = 0;
    for i = 1:input
        for j = 1:input
        state = state + i*1e-6;
        end
    end
end