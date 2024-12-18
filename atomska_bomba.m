function [t_1_opt, fi_1_opt, t_3_opt, d_opt] = atomska_bomba(primer)

%citanje ulaznih podataka
run(primer);

%racunanje tilde_t
tilde_t = r*pi./v_a;

%veza ugla i t_1
fi_1 = @(t_1) v_a*t_1./r;

%polozaj bombe u trenutku t_1
x_1 = @(t_1) r.*sin(fi_1(t_1));
y_1 = @(t_1) r.*(1-cos(fi_1(t_1)));

%brzina bombe u trenutku t_1
v_1x = @(t_1) v_a.*cos(fi_1(t_1));
v_1y = @(t_1) v_a.*sin(fi_1(t_1));

%trenutak u kom bomba eksplodira
t_2 = @(t_1) t_1+(v_1y(t_1)+sqrt(v_1y(t_1).^2+2*g.*(y_1(t_1)+h)))./g;

%polozaj bombe u trenutku t_2
x2 = @(t_1) x_1(t_1)+v_1x(t_1).*(t_2(t_1)-t_1);

%koeficijenti u kvadratnoj jednacini po t_3
k_2 = v_t.^2-v_a.^2;
k_1 = @(t_1) 2*v_a.*(x2(t_1)+v_a.*(t_2(t_1)-tilde_t));
k_0 = @(t_1) (x2(t_1)+v_a.*(t_2(t_1)-tilde_t)).^2+(h+2*r).^2;

%racunanje t_3
t_3 = @(t_1) (k_1(t_1)+sqrt(k_1(t_1).^2+4*k_2.*k_0(t_1)))./(2*k_2);

%priprema za iscrtavanje grafika
X = linspace(0,tilde_t);
figure;
screen_size = get(0, 'ScreenSize');
set(gcf, 'Position', [(screen_size(3)-1080)/2, (screen_size(4)-540)/2, 1080, 540]);

%prvi grafik
subplot(1,2,1);
plot(X, t_3(X));
set(gca,'position',[0.05, 0.1, 0.425, 0.825]);
ylim([0, 300]);
xlim([0, tilde_t]); 
title('Зависност t_{3} од t_{1}');
xlabel('$t_{1}$[s]', 'Interpreter', 'latex');
ylabel('$t_{3}$[s]', 'Interpreter', 'latex');

%drugi grafik
subplot(1,2,2)
plot(fi_1(X), t_3(X));
set(gca,'position',[0.55, 0.1, 0.425, 0.825]);
ylim([0, 300]);
xlim([0, pi]); 
title('Зависност t_{3} од φ_{1}');
xlabel('$\varphi_{1}$[rad]', 'Interpreter', 'latex');
ylabel('$t_{3}[s]$', 'Interpreter', 'latex');

%vracanje optimalnih vrednosti
t_1_opt = fminbnd(@(t_1) -t_3(t_1),0,tilde_t);
fi_1_opt = fi_1(t_1_opt);
t_3_opt = t_3(t_1_opt);
d_opt = v_t.*t_3_opt;
