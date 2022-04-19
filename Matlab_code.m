%% EXAMPLE NUMERIC SOLVER FOR EQUILIBRIUM %%
syms x
a = 10;
S = vpasolve(a*x+3==0,x)

eq1 = a*x+3==0;
a = 3;
S = vpasolve(eq1,x)


eq1 = a*x+3==0;
S = vpasolve(eq1,x)
% What this teaches us is that we need to re-read in the equation every
% time we change a parameter value


%% Basic solution and learning how to manipulate it

syms CL CD MS MC

gC = 0.1;
dC = 0.001;
gammaC = 0.001;
gM = 0.6;
gMCL = 0.06;
dM = 0.05;
hMS = 0.3;
gammaM = 0.001;
er = 0.01;
ref = 0.5;

eq1 = (gC * CL + gammaC) * (1 - CL - CD - MS - MC) - gMCL * CL * (MS + MC) - dC * CL==0;
eq2 = dC * CL + MC * (hMS * (1 - ref) + dM) - gM * CD * (MS + MC) - er * CD -  gammaM * CD==0;
eq3 = (gM * MS + gM * MC + gammaM) * (1 - CL - CD - MS - MC) - MS * (hMS + dM)==0;
eq4 = gM * (MC + MS) * CD + gMCL * CL * (MS + MC) + gammaM * CD - MC * (hMS * (1 - ref) + dM + er)==0;

S = vpasolve([eq1, eq2, eq3, eq4],[CL, CD, MS, MC]); % VPASOLVE finds equilibria numerically
S.CL % Shows, as a vertical vector, the equilibrium values of live coral that are possible

res = [S.CL,S.CD,S.MS,S.MC]'; % stores the equilibria in a matrix. We use the transpose here so that the rows are state variables, and the columns are equilibria
min(res) % computes minima in each column
min(res)>=0 % finds which of these minima are non-negative; this allows us to focus on non-negative equilibria only
res(:,min(res)>=0) % extracts just these equilibria
sort(res(1,min(res)>=0)) % extracts just the live coral abundances for these equilibria and sorts them from smallest to largest


%% Learning more about how MATLAB plots surfaces

figure(1)
xset = [0,1,2,3,4,0];
yset = [0,1,2,3,4,0];
zset = [0,1,2,3,4,5];
plot3(xset,yset,zset,'linestyle','none','marker','o') 
% OK, so one of the first things we can do is just concatenate all this
% stuff together and make a 3D plot of points

figure(2)
xset = [0,1];
yset = [0,1];
zset = [0,0;0,0];
surf(xset,yset,zset)
hold on
surf(xset,yset,zset+1)

figure(3)
plot3(xset,yset,zset,'linestyle','none','marker','o')

% Alternatively, we can put the high, low, and unstable equilibria on
% different surfaces and plot those layered (using hold on)

%% Let's GOOOOOOO 

refset = linspace(0,1,200);
herbset = linspace(0,1,205);

highstab = NaN(length(refset),length(herbset)); % Placeholder for high stabile equilibrium
unstab = highstab; % Placeholder for unstable equilibrium
lowstab = highstab; % Placeholder for low unstable equilibrium

%%

% to run still:
% i = 225, j > 6
% i = 248, j > 52
% i = 265, j > 283
tic
for i = 249:length(herbset)
    hMS = herbset(i) % update herbivory value
    eqctr = 0; % Re-set the equilibrium counter
    for j = 1:length(refset)
        ref = refset(j); % update refuge value
        eq1 = (gC * CL + gammaC) * (1 - CL - CD - MS - MC) - gMCL * CL * (MS + MC) - dC * CL==0;
        eq2 = dC * CL + MC * (hMS * (1 - ref) + dM) - gM * CD * (MS + MC) - er * CD -  gammaM * CD==0;
        eq3 = (gM * MS + gM * MC + gammaM) * (1 - CL - CD - MS - MC) - MS * (hMS + dM)==0;
        eq4 = gM * (MC + MS) * CD + gMCL * CL * (MS + MC) + gammaM * CD - MC * (hMS * (1 - ref) + dM + er)==0;

        S = vpasolve([eq1, eq2, eq3, eq4],[CL, CD, MS, MC]); % VPASOLVE finds equilibria numerically
        res = [S.CL,S.CD,S.MS,S.MC]'; % Store the output as a matrix
        res = res(:,abs(sum(imag(res)))==0); % Remove all imaginary roots
        CLstar = sort(res(1,min(res)>=0)); % Remove all negative roots
        
        if length(CLstar)>1 % If there are multiple equilibria, we fill them from high to low. 
            % Note that this represents an ASSUMPTION that the MIDDLE
            % equilibrium is UNSTABLE. Based on the Mathematica script
            % (which includes jacobian matrix analysis), this is true for
            % LIVE CORAL. However, it is not always true for other state
            % variables (e.g., dead coral), so use caution!
            if length(CLstar) == 2 % if there are exactly 2 equilibria one's going to be semi-stable (on the boundary w/ the unstable equilibrium)
                highstab(j,i) = CLstar(2);
                lowstab(j,i) = CLstar(1);
                %if eqctr == 0 % if we're just starting to enter the region of bistability
                %    unstab(j,i) = CLstar(2); % then the unstable equilibrium is the same as the high-coral stable one
                %else % Otherwise, it's the same as the low-coral stable one
                %    unstab(j,i) = CLstar(1);
                %end
            else
                eqctr = 1; % Once we get to the 3-equilibrium state, set eqctr to 1
                highstab(j,i) = CLstar(3);
                unstab(j,i) = CLstar(2);
                lowstab(j,i) = CLstar(1);
            end
        else % If there's only one equilibrium, we have to decide if it's the 'high' one or the 'low' one
            if eqctr == 0
                % Remember how we set eqctr == 0 outside of this for loop?
                % Within THIS for loop (the one indexed by j), we're
                % iterating over a range of refuge values. When the algal
                % refuge = 0, we expect to be at the 'high' coral
                % equilibrium (or, at least, the highest value we can be
                % at). So we'll save the sole equilibrium here.
                
                % However, when there's absolutely no herbivory, sometimes
                % we never find ourselves at the high equilibrium
                % initially. So we'll put in a "soft check" here that if
                % the equilibrium coral cover is less than 0.2, you call it
                % low-coral anyway
                if CLstar > 0.2
                    highstab(j,i) = CLstar;
                else
                    lowstab(j,i) = CLstar;
                end
            else
                % Otherwise, when we're 'past' the region of bistability
                % (that is, eqctr has now been set to 1 because we've gone
                % through a region with more than 1 equilibrium), we're
                % going to save the equilibrium as the "low coral" state.
                lowstab(j,i) = CLstar;
            end
        
        end
        
    end
    toc
    save('500by505run')
end
%%
%%%save('50by70run')
%%%save('200by205run')
%%

figure(1)
clf(1)
surf(herbset,refset,highstab)
hold on
surf(herbset,refset,lowstab)
surf(herbset,refset,unstab)
xlabel('Strength of Herbivory')
ylabel('Magnitude of Algal Refuge')
zlabel('Coral Cover')
meshgrid off
shading interp
%view(70,35)

%% How do we smooth the edges?
% We could add a little extra around the edges of the unstable equilibrium,
% creating overlap with the stable ones on either side

unstab_mod = unstab; % Create a new placeholder matrix for unstable equilibria
lowstab_mod = lowstab;
highstab_mod = highstab;

for i = 1:length(herbset)
    eqctr = 0; % start the equilibrium counter at 0
    for j = 1:length(refset)
        
        % Check the "high equilibrium" boundary
        if j < length(refset) - 1 % be careful not to go off the edge!
            if isnan(unstab(j,i))
                if isnan(unstab(j+1,i))
                else % if the next-in-line unstable equilibrium exists
                    unstab_mod(j+1,i) = lowstab(j+1,i);
                    unstab_mod(j+2,i) = lowstab(j+2,i);
                end
            else
                eqctr = 1; % once you enter the region where you have unstable equilibria, turn this on
            end
        end
        
        % Check the "low equilibrium" boundary
        if j > 2 % be careful not to go off the edge!
            if isnan(unstab(j,i))
                if isnan(unstab(j-1,i))
                else % If the just-previous unstable equilibrium exists
                    unstab_mod(j-1,i) = highstab(j-1,i);
                    unstab_mod(j-2,i) = highstab(j-2,i);
                end
            end
        end
        
    end
end

%% How do we smooth the edges?
% We could add a little extra around the edges of the unstable equilibrium,
% creating overlap with the stable ones on either side

unstab_mod = unstab; % Create a new placeholder matrix for unstable equilibria
lowstab_mod = lowstab;
highstab_mod = highstab;

for i = 1:length(herbset)
    eqctr = 0; % start the equilibrium counter at 0
    for j = 1:length(refset)
        
        % Check the "high equilibrium" boundary
        if j < length(refset) - 1 % be careful not to go off the edge!
            if isnan(unstab(j,i))
                if isnan(unstab(j+1,i))
                else % if the next-in-line unstable equilibrium exists
                    unstab_mod(j+1,i) = lowstab(j+1,i);
                    unstab_mod(j+2,i) = lowstab(j+2,i);
                    unstab_mod(j+3,i) = lowstab(j+3,i);
                        unstab_mod(j+4,i) = lowstab(j+4,i);
                    if j < length(refset)/2
                        unstab_mod(j+5,i) = lowstab(j+5,i);
                        unstab_mod(j+6,i) = lowstab(j+6,i);
                        if j < length(refset)/5
                        unstab_mod(j+7,i) = lowstab(j+7,i);
                        unstab_mod(j+8,i) = lowstab(j+8,i);
                        unstab_mod(j+9,i) = lowstab(j+9,i);
                        end
                            
                    end
                end
            else
                eqctr = 1; % once you enter the region where you have unstable equilibria, turn this on
            end
        end
        
        % Check the "low equilibrium" boundary
        if j > 2 % be careful not to go off the edge!
            if isnan(unstab(j,i))
                if isnan(unstab(j-1,i))
                else % If the just-previous unstable equilibrium exists
                    unstab_mod(j-1,i) = highstab(j-1,i);
                    unstab_mod(j-2,i) = highstab(j-2,i);
                    unstab_mod(j-3,i) = highstab(j-3,i);
                        unstab_mod(j-4,i) = highstab(j-4,i);
                    if j < length(refset)/2
                        unstab_mod(j-5,i) = highstab(j-5,i);
                        unstab_mod(j-6,i) = highstab(j-6,i);
                        unstab_mod(j-7,i) = highstab(j-7,i);
                        unstab_mod(j-8,i) = highstab(j-8,i);
                    end
                end
            end
        end
        
    end
end

%%
figure(2)
clf(2)
surf(herbset,refset,highstab,'EdgeColor','n','FaceLighting','flat')
hold on
surf(herbset,refset,lowstab)
surf(herbset,refset,unstab_mod,'FaceColor','red','EdgeColor','n')
xlabel('Strength of Herbivory')
ylabel('Magnitude of Algal Refuge')
zlabel('Coral Cover')
meshgrid off
shading interp
view(-165,25)
lighting none

% plot no-hysteresis line
nohystchoice = 45;
plot3(herbset,refset(nohystchoice)*ones(length(herbset)),lowstab(nohystchoice,:))


%%
figure(2)
close(2)
figure(2)
coords = get(gcf,'Position');
set(gcf,'Position',[coords(1),coords(2),coords(3)*1.5,coords(4)])

herblim = 0.5;

subplot(2,7,[1:4,8:11])
surf(refset,herbset,highstab','EdgeColor','n','FaceLighting','flat')
hold on
surf(refset,herbset,lowstab')
surf(refset,herbset,unstab_mod','FaceColor','red','EdgeColor','n')
ylabel('Herbivores')
xlabel('Algal Refugia')
zlabel('Coral Cover')
meshgrid off
shading interp
view(70,35)
ylim([0,herblim])
%lighting none

% plot no-hysteresis line
nohystchoice = round(length(refset)*3/4,0);
nohystchoice = round(length(refset)*9/10,0);
plot3(refset(nohystchoice)*ones(length(herbset)),herbset,lowstab(nohystchoice,:),'Color','k','LineWidth',3)

% plot hysteresis line
hystchoice = round(length(refset)/5,0);

% select lower edge
lowhyst = lowstab(hystchoice,:);
lowhyst = lowhyst(isfinite(lowhyst));
lowhystx = herbset(isfinite(lowstab(hystchoice,:)));

% select unstable edge
unhyst = unstab(hystchoice,:);
unhyst = unhyst(isfinite(unhyst));
unhystx = herbset(isfinite(unstab(hystchoice,:)));

% select upper edge
highhyst = highstab(hystchoice,:);
highhyst = highhyst(isfinite(highhyst));
highhystx = herbset(isfinite(highstab(hystchoice,:)));

% concatenate and plot
yset = [lowhystx,wrev(unhystx),highhystx];
xset = refset(hystchoice)*ones(length(yset));
zset = [lowhyst,wrev(unhyst),highhyst];
plot3(xset,yset,zset,'LineWidth',3,'Color',[0.6350 0.0780 0.1840])


% ADD ARROW
shiftchoice = round(length(herbset)/10,0);
plot3(refset(hystchoice:nohystchoice),herbset(shiftchoice)*ones(nohystchoice-hystchoice+1),lowstab(hystchoice:nohystchoice,shiftchoice),'LineStyle',':','LineWidth',3,'Color','white')
plot3([refset(nohystchoice)-.02,refset(nohystchoice)],[herbset(shiftchoice)-.02,herbset(shiftchoice)],[lowstab(nohystchoice)+.05,lowstab(nohystchoice)+.01],'LineWidth',3,'Color','white')
plot3([refset(nohystchoice)-.04,refset(nohystchoice)],[herbset(shiftchoice)+.005,herbset(shiftchoice)],[lowstab(nohystchoice)+.05,lowstab(nohystchoice)+.01],'LineWidth',3,'Color','white')

%plot3(refset(nohystchoice),herbset(shiftchoice),lowstab(nohystchoice,shiftchoice),'>','LineWidth',3,'Color','white')

text(0.05,0.01,.95,'A','FontSize',20)
text(refset(nohystchoice),herbset(shiftchoice)+.01,lowstab(nohystchoice,shiftchoice)+.05,'C','FontSize',20,'Color','white')
text(refset(hystchoice),herbset(shiftchoice)-.03,lowstab(nohystchoice,shiftchoice)+.05,'B','FontSize',20,'Color','white')

subplot(2,6,[5:6])
%plot(yset,zset,'LineWidth',3,'Color',[0.6350 0.0780 0.1840])
hold on
plot(lowhystx,lowhyst,'LineWidth',3,'Color',[0.6350 0.0780 0.1840])
plot(highhystx,highhyst,'LineWidth',3,'Color',[0.6350 0.0780 0.1840])
plot(unhystx,unhyst,'LineWidth',3,'Color',[0.6350 0.0780 0.1840],'LineStyle',':')
xlabel('Herbivores')
ylabel('Coral Cover')
ylim([0,1])
text(0.01,0.91,'B','FontSize',20)
xlim([0,herblim])

subplot(2,6,[11:12])
plot(herbset,lowstab(nohystchoice,:),'LineWidth',3,'Color','k')
xlabel('Herbivores')
ylabel('Coral Cover')
ylim([0,1])
text(0.01,0.91,'C','FontSize',20)
xlim([0,herblim])

%% Try next

% Can we plot surf on a list of points instead of a matrix

% What about contour plots

% What about turning the system "on its side" and solving for the refuge
% level given C_L?


% Change the colours of the surface and add the lines representing the
% high-refuge and low-refuge scenarios


%% Updated figure -- creating bifurcation data
refset2 = linspace(0,1,100);

hystchoice = round(length(herbset)/10*3,0);
hMS = herbset(hystchoice); % update herbivory value

CLstar_high = NaN*refset2;
CLstar_low = CLstar_high;
CLstar_un = CLstar_high;

Mstar_high = NaN*refset2;
Mstar_low = NaN*refset2;
Mstar_un = NaN*refset2;

eqctr = 0;
    for j = 1:length(refset2)
        ref = refset2(j); % update refuge value
        eq1 = (gC * CL + gammaC) * (1 - CL - CD - MS - MC) - gMCL * CL * (MS + MC) - dC * CL==0;
        eq2 = dC * CL + MC * (hMS * (1 - ref) + dM) - gM * CD * (MS + MC) - er * CD -  gammaM * CD==0;
        eq3 = (gM * MS + gM * MC + gammaM) * (1 - CL - CD - MS - MC) - MS * (hMS + dM)==0;
        eq4 = gM * (MC + MS) * CD + gMCL * CL * (MS + MC) + gammaM * CD - MC * (hMS * (1 - ref) + dM + er)==0;

        S = vpasolve([eq1, eq2, eq3, eq4],[CL, CD, MS, MC]); % VPASOLVE finds equilibria numerically
        res = [S.CL,S.CD,S.MS,S.MC]'; % Store the output as a matrix
        res = res(:,abs(sum(imag(res)))==0); % Remove all imaginary roots
        res = res(:,min(res)>=0); % extracts just these equilibria

if size(res,2) == 1
    if eqctr == 1
        CLstar_low(j) = res(1,1);
        Mstar_low(j) = res(3,1)+res(4,1);
    else
        CLstar_high(j) = res(1,1);
        Mstar_high(j) = res(3,1)+ res(4,1);
    end
else
    eqctr = 1;
    CLstars = sort(res(1,:));
    Mstars = res(3,:)+res(4,:);
    Mstars = sort(Mstars);
    CLstar_low(j) = CLstars(1);
    CLstar_un(j) = CLstars(2);
    CLstar_high(j) = CLstars(3);
    Mstar_low(j) = Mstars(3);
    Mstar_un(j) = Mstars(2);
    Mstar_high(j) = Mstars(1);
    
end
        
    end
    
%% Make the updated plot
figure(2)
close(2)
figure(2)
coords = get(gcf,'Position');
set(gcf,'Position',[coords(1),coords(2),coords(3)*1.5,coords(4)])

herblim = 0.4;
corcol = [107/255,126/255,158/255];
maccol = [144/255,126/255,81/255];

subplot(2,6,[1:2])
plot(refset2,CLstar_high,'LineWidth',2,'Color',corcol)
hold on
plot(refset2,CLstar_low,'LineWidth',2,'Color',corcol)
plot(refset2,CLstar_un,'LineStyle','--','LineWidth',2,'Color',corcol)
xlabel('Algal Protection (\rho)')
ylabel('Coral (C_L)')
ylim([0,1])

subplot(2,6,[7:8])
plot(refset2,Mstar_high,'LineWidth',2,'Color',maccol)
hold on
plot(refset2,Mstar_low,'LineWidth',2,'Color',maccol)
plot(refset2,Mstar_un,'LineStyle','--','LineWidth',2,'Color',maccol)
xlabel('Algal Protection (\rho)')
ylabel('Macroalgae (M_S+M_C)')
ylim([0,1])

subplot(2,7,[4:7,11:14])
surf(refset,herbset,highstab','EdgeColor','n','FaceLighting','flat')
hold on
surf(refset,herbset,lowstab')
surf(refset,herbset,unstab_mod','FaceColor','red','EdgeColor','n')
ylabel('Herbivores (h)')
xlabel('Algal Protection (\rho)')
zlabel('Coral Cover (C_L)')
meshgrid off
shading interp
%view(70,35)
view(33,30)
ylim([0,herblim])
%lighting none
colormap(flipud(parula))

% plot hysteresis line
hystchoice = round(length(herbset)/10*3,0);

% select lower edge
lowhyst = lowstab(:,hystchoice);
lowhyst = lowhyst(isfinite(lowhyst));
lowhystx = refset(isfinite(lowstab(:,hystchoice)));

% select unstable edge
unhyst = unstab(:,hystchoice);
unhyst = unhyst(isfinite(unhyst));
unhystx = refset(isfinite(unstab(:,hystchoice)));

% select upper edge
highhyst = highstab(:,hystchoice);
highhyst = highhyst(isfinite(highhyst));
highhystx = refset(isfinite(highstab(:,hystchoice)));

% concatenate and plot
xset = [highhystx,wrev(unhystx),lowhystx];
yset = herbset(hystchoice)*ones(length(xset));
zset = [highhyst',wrev(unhyst'),lowhyst'];
%plot3(xset,yset,zset,'LineWidth',3,'Color',[0.6350 0.0780 0.1840])
plot3(xset,yset,zset,'LineWidth',3,'Color','k')






