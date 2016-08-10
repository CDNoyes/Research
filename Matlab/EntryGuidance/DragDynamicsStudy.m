function DragDynamicsStudy()

load('EntryGuidance/HighElevationPlanner/Trajectories/ReferenceTrajectory_DR780_CR0.mat');


delta_CD = linspace(-.4,.4,16);
delta_rho_per = linspace(-.3,.3,17);
for k = 1:length(delta_CD)
    for j = 1:length(delta_rho_per)
        for i = 1:length(ref.time)
            delta_rho(j) = delta_rho_per(j)*ref.rho(i);
            D = ref.D(i) + delta_CD(k)*ref.D(i)/ref.C_D(i) + ref.D(i)*delta_rho(j)/ref.rho(i) + ref.D(i)/ref.C_D(i)*delta_rho(j)/ref.rho(i)*delta_CD(k);
            [a(i),b(i)] = DragFBL(ref.g(i),ref.L(i),D,ref.state(i,1),ref.state(i,4),ref.state(i,5),ref.rho(i)+delta_rho(j),ref.rho_dot(i),[],ref.C_D(i)+delta_CD(k),ref.C_D_dot(i));
            a_max(k,j) = max((a));
            b_max(k,j) = max((b));
            a_min(k,j) = min(a);
            b_min(k,j) = min(b);
        end
    end
end

figure
surf(delta_CD,delta_rho_per,a_max')
hold on
surf(delta_CD,delta_rho_per,a_min')
xlabel('C_D offset (-)')
ylabel('\rho offset(%)')
zlabel('a')
title('Bounds on drag dynamic term a')

figure
surf(delta_CD,delta_rho_per,b_max')
hold on
surf(delta_CD,delta_rho_per,b_min')
xlabel('C_D offset (-)')
ylabel('\rho offset(%)')
zlabel('b')
title('Bounds on drag dynamic term b')

end