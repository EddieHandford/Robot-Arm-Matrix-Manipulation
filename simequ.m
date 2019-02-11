
z = 1.4650;
y = 1.3205;

syms za

eqn1 = -sqrt(1-(z-1.60131977*sin(za+0.040603))/0.6)+2.668866167*cos(za-0.040603) == y;


sol = solve([eqn1] , [za]);
zaSol= sol.za
