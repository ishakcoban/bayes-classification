load('classes.mat')

plot(C1(1,:), C1(2,:),'+')
hold
plot(C2(1,:), C2(2,:),'rd')
plot(C3(1,:), C3(2,:),'*')

checkClassSizes(C1,C2,C3);

classes = {'C1','C2','C3'};
classMatrixes = {C1,C2,C3};

x1 = [6,8]';
x2 = [-5,-3]';
x3 = [12,0]';

vectors = {x1,x2,x3};

for i = 1:length(vectors)

    classValues = {0,0,0};
    for j = 1:length(classes)
        m = mean(transpose(classMatrixes{1,j}));
        S = cov(classMatrixes{1,j}');
        x = vectors{1,i};
        classValues{1,j} =  comp_gauss_dens_val(m',S,x);
    end

    disp("x" + i +"_pdf = "+ classValues{1,1} + "_for_C1, "+ classValues{1,2} + "_for_C2, "+ classValues{1,3} + "_for_C3");

    %Q5b
    maxValue = max(cell2mat(classValues));
    for j = 1:length(classValues)
        if(maxValue == classValues{1,j})
            disp("x" + i + " belongs to class " + classes{1,j});
        end
    end
    disp(" ");
end

function checkClassSizes(in1, in2, in3)

    if size(in1) == size(in2) & size(in2) == size(in3) & size(in1) == size(in3)
        disp('Classes are equal');
    else
        disp('Classes are not equal');
    end

end

function z = comp_gauss_dens_val(m,S,x)

    [l,q]=size(m);
	    z= (1/((2*pi)^ (l/2)*det(S)^ 0.5)) ...
		    *exp(-0.5*(x-m)'*inv(S)*(x-m));

end
