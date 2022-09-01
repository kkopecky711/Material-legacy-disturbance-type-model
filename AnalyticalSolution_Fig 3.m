%% Basic solution and learning how to manipulate it

% This code shows how you can use a semi-analytical approach in which you
% ask MATLAB to generate the numerical values of the equilibria for a set
% of parameter values.


% First, define the symbols
syms CL CP CD ME MP

% Second, define the parameters
gC = 0.1;
dC = 0.001;
gammaC = 0.001;
s = 1;
gM = 0.6;
gMC = 0.06;
dM = 0.05;
h = 0.3;
gammaM = 0.001;
er = 0.01;
p = 0.5;

% Third, input the system of five differential equations. Note that here
% we've made the substitution that S = 1 - CL - CP - CD - ME - MP.
eq1 = (gC * CL + gammaC) * (1 - CL - CP - CD - ME - MP) - gMC * CL * (ME + MP) - dC * CL==0;
eq2 = gC*CP*CD + s*gammaC*CD - gMC*CP*(ME+MP) - dC*CP - er*CP == 0;
eq3 = dC * (CL + CP) + MP * (h * (1 - p) + dM) - gC*CP*CD -  s*gammaC*CD - gM * CD * (ME + MP) - er * CD -  gammaM * CD==0;
eq4 = (gM * ME + gM * MP + gammaM) * (1 - CL - CP - CD - ME - MP) - ME * (h + dM)==0;
eq5 = gM * (ME + MP) * CD + gMC * (CL+CP) * (ME + MP) + gammaM * CD - MP * (h * (1 - p) + dM + er)==0;

% Use a numerical solver to find the exact values of the equilibria
S = vpasolve([eq1, eq2, eq3, eq4, eq5],[CL, CP, CD, ME, MP]); % VPASOLVE finds equilibria numerically
S.CL % Shows, as a vertical vector, the equilibrium values of live coral that are possible

% Extract the relevant results; here are just some examples of how you can
% examine them. Note that we are only interested in real, positive
% equilibria, so you'll see that we subset for these below.
res = [S.CL,S.CP,S.CD,S.ME,S.MP]'; % stores the equilibria in a matrix. We use the transpose here so that the rows are state variables, and the columns are equilibria
res = res(:,abs(sum(imag(res)))==0); % removes imaginary equilibria
min(res) % computes minima in each column
min(res)>=0 % finds which of these minima are non-negative; this allows us to focus on non-negative equilibria only
res(:,min(res)>=0) % extracts just these equilibria
sort(res(1,min(res)>=0)) % extracts just the live coral abundances for these equilibria and sorts them from smallest to largest


%% Large simulation to generate surface
% Vectors of parameter values to iterate over
refset = linspace(0,1,10);
herbset = linspace(0,1,15);

% Set up efficient data collection by creating placeholder matrices
highstab_CL = NaN(length(refset),length(herbset)); % Placeholder for high stable equilibrium
unstab_CL = highstab_CL; % Placeholder for unstable equilibrium
lowstab_CL = highstab_CL; % Placeholder for low unstable equilibrium

% placeholders for macroalgae; note that we have to be careful to sum the
% macroalgal abundances before saving them.
highstab_M = highstab_CL;
unstab_M = highstab_CL;
lowstab_M = highstab_CL;

%%

% Collect data on equilibria
tic % sets a timer; useful to understand the performance of this loop.
for i = 1:length(herbset) % Loop over every herbivory value
    h = herbset(i) % update herbivory value
    eqctr = 0; % Re-set the equilibrium counter
    for j = 1:length(refset) % Loop over every refuge strength
        p = refset(j); % update refuge value
        
        % Use the routine explored above to find the real, nonnegative
        % equilibria
        eq1 = (gC * CL + gammaC) * (1 - CL - CP - CD - ME - MP) - gMC * CL * (ME + MP) - dC * CL==0;
        eq2 = gC*CP*CD + s*gammaC*CD - gMC*CP*(ME+MP) - dC*CP - er*CP == 0;
        eq3 = dC * (CL + CP) + MP * (h * (1 - p) + dM) -  s*gammaC*CD - gC*CP*CD - gM * CD * (ME + MP) - er * CD -  gammaM * CD==0;
        eq4 = (gM * ME + gM * MP + gammaM) * (1 - CL - CP - CD - ME - MP) - ME * (h + dM)==0;
        eq5 = gM * (ME + MP) * CD + gMC * (CL+CP) * (ME + MP) + gammaM * CD - MP * (h * (1 - p) + dM + er)==0;

        S = vpasolve([eq1, eq2, eq3, eq4, eq5],[CL, CP, CD, ME, MP]); % VPASOLVE finds equilibria numerically
        res = [S.CL,S.CP,S.CD,S.ME,S.MP]'; % stores the equilibria in a matrix. We use the transpose here so that the rows are state variables, and the columns are equilibria
        res = res(:,abs(sum(imag(res)))==0); % Remove all imaginary roots
        res = res(:,min(res)>=0); % Remove all negative roots
        CLstar = sort(res(1,:)+res(2,:)); % Sort and save the Coral values
        Mstar = sort(res(4,:)+res(5,:)); % Sort and save the Macroalgae values
        
        % Save the equilibria in their appropriate holding matrices
        if length(CLstar)>1 % If there are multiple equilibria, we fill them from high to low. 
            % Note that this represents an ASSUMPTION that the MIDDLE
            % equilibrium is UNSTABLE. Based on the Mathematica script
            % (which includes jacobian matrix analysis), this is true for
            % LIVE CORAL. However, it is not always true for other state
            % variables (e.g., dead coral), so use caution!
            if length(CLstar) == 2 % if there are exactly 2 equilibria one's going to be semi-stable (on the boundary w/ the unstable equilibrium)
                highstab_CL(j,i) = CLstar(2);
                lowstab_CL(j,i) = CLstar(1);
                %if eqctr == 0 % if we're just starting to enter the region of bistability
                %    unstab(j,i) = CLstar(2); % then the unstable equilibrium is the same as the high-coral stable one
                %else % Otherwise, it's the same as the low-coral stable one
                %    unstab(j,i) = CLstar(1);
                %end
                
                lowstab_M(j,i) = Mstar(1);
                highstab_M(j,i) = Mstar(2);
                
            else
                eqctr = 1; % Once we get to the 3-equilibrium state, set eqctr to 1
                highstab_CL(j,i) = CLstar(3);
                unstab_CL(j,i) = CLstar(2);
                lowstab_CL(j,i) = CLstar(1);
                
                highstab_M(j,i) = Mstar(3);
                unstab_M(j,i) = Mstar(2);
                lowstab_M(j,i) = Mstar(1);
                
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
                    highstab_CL(j,i) = CLstar;
                    lowstab_M(j,i) = Mstar;
                else
                    lowstab_CL(j,i) = CLstar;
                    highstab_M(j,i) = Mstar;
                end
            else
                % Otherwise, when we're 'past' the region of bistability
                % (that is, eqctr has now been set to 1 because we've gone
                % through a region with more than 1 equilibrium), we're
                % going to save the equilibrium as the "low coral" state.
                lowstab_CL(j,i) = CLstar;
                highstab_M(j,i) = Mstar;
            end
        
        end
        
    end
    toc % Reports out the time that iterating through this value went
    save('test') 
    % ^Highly recommend saving at each iteration of the loop, as this is a 
    % computationally intensive analysis, and it'd be a shame to lose the
    % data. But note that you want to be very careful to not over-write
    % your matrices with zeros, which is why this piece of code is in a
    % separate chunk.
end

%% Modify data to address numerical errors

% Simulation outputs that were erroneous were very rare (one in 10000) and
% extremely obvious (because they predicted coral equilibria >> 1. We
% simply replaced these with NaNs to prevent them from impacting the
% plotting of the results.
for i = 1:length(herbset)
    for j = 1:length(refset)
        if highstab_CL(j,i) > 1
            highstab_CL(j,i) = NaN;
            highstab_M(j,i) = NaN;
        end
    end
end

% Remove missing data, if any.
highstab_CL(:,103) = [];
lowstab_CL(:,103) = [];
unstab_CL(:,103) = [];
highstab_M(:,103) = [];
lowstab_M(:,103) = [];
unstab_M(:,103) = [];
herbset(103) = [];


%% Smoothing edges for plotting a 3D surface
% A limitation of MATLAB is that its "surf" function cannot plot a surface
% that folds back on itself, as this landscape does. Therefore, if we plot
% the unmodified surfaces, we wind up with little gaps between the two
% stable equilibria and the unstable equilibrium. To fix this, we "patch"
% the diagram by adding a little extra around the edges of the unstable
% equilibrium, creating overlap with the stable ones on either side

unstab_mod_CL = unstab_CL; % Create a new placeholder matrix for unstable equilibria
unstab_mod_M = unstab_M; % Create a new placeholder matrix for unstable equilibria

for i = 1:length(herbset)
    eqctr = 0; % start the equilibrium counter at 0
    for j = 1:length(refset)
        
        % Check the "high equilibrium" boundary
        if j < length(refset) - 1 % be careful not to go off the edge!
            if isnan(unstab_CL(j,i))
                if isnan(unstab_CL(j+1,i))
                else % if the next-in-line unstable equilibrium exists
                    unstab_mod_CL(j+1,i) = lowstab_CL(j+1,i);
                    unstab_mod_CL(j+2,i) = lowstab_CL(j+2,i);
                    unstab_mod_CL(j+3,i) = lowstab_CL(j+3,i);
                        unstab_mod_CL(j+4,i) = lowstab_CL(j+4,i);
                    if j < length(refset)/2
                        unstab_mod_CL(j+5,i) = lowstab_CL(j+5,i);
                        unstab_mod_CL(j+6,i) = lowstab_CL(j+6,i);
                        if j < length(refset)/3
                        unstab_mod_CL(j+7,i) = lowstab_CL(j+7,i);
                        unstab_mod_CL(j+8,i) = lowstab_CL(j+8,i);
                        unstab_mod_CL(j+9,i) = lowstab_CL(j+9,i);
                        unstab_mod_CL(j+10,i) = lowstab_CL(j+10,i);
                        unstab_mod_CL(j+11,i) = lowstab_CL(j+11,i);
                            if j < length(refset)/5
                                unstab_mod_CL(j+12,i) = lowstab_CL(j+12,i);
                                unstab_mod_CL(j+13,i) = lowstab_CL(j+13,i);
                            end
                        end
                            
                    end
                end
            else
                eqctr = 1; % once you enter the region where you have unstable equilibria, turn this on
            end
        end
        
        % Check the "low equilibrium" boundary
        if j > 2 % be careful not to go off the edge!
            if isnan(unstab_CL(j,i))
                if isnan(unstab_CL(j-1,i))
                else % If the just-previous unstable equilibrium exists
                    unstab_mod_CL(j-1,i) = highstab_CL(j-1,i);
                    unstab_mod_CL(j-2,i) = highstab_CL(j-2,i);
                    unstab_mod_CL(j-3,i) = highstab_CL(j-3,i);
                        unstab_mod_CL(j-4,i) = highstab_CL(j-4,i);
                    if j < length(refset)/2
                        unstab_mod_CL(j-5,i) = highstab_CL(j-5,i);
                        unstab_mod_CL(j-6,i) = highstab_CL(j-6,i);
                        unstab_mod_CL(j-7,i) = highstab_CL(j-7,i);
                        unstab_mod_CL(j-8,i) = highstab_CL(j-8,i);
                        if j < length(refset)/5
                        unstab_mod_CL(j-9,i) = highstab_CL(j-9,i);
                        unstab_mod_CL(j-10,i) = highstab_CL(j-10,i);
                        unstab_mod_CL(j-11,i) = highstab_CL(j-11,i);
                        %unstab_mod_CL(j-12,i) = highstab_CL(j-12,i);
                        end
                    end
                end
            end
        end
        
    end
end

%% Plot the surfaces

figure(2)
close(2)
figure(2)
coords = get(gcf,'Position');
set(gcf,'Position',[coords(1),coords(2),coords(3)*2,coords(4)*1.1])

subplot(1,2,1)
surf(herbset,refset,highstab_CL,'EdgeColor','n','FaceLighting','flat')
hold on
surf(herbset,refset,lowstab_CL)
surf(herbset,refset,unstab_mod_CL,'FaceColor','red','EdgeColor','n')
xlabel('Strength of Herbivory')
ylabel('Magnitude of Algal Refuge')
zlabel('Coral Cover')
meshgrid off
shading interp
view(-53,28)
%view(33,30)
%view(70,35)
lighting none


subplot(1,2,2)
surf(herbset,refset,highstab_M,'EdgeColor','n','FaceLighting','flat')
hold on
surf(herbset,refset,lowstab_M)
surf(herbset,refset,unstab_M,'FaceColor','red','EdgeColor','n')
xlabel('Strength of Herbivory')
ylabel('Magnitude of Algal Refuge')
zlabel('Macroalgal Cover')
meshgrid off
shading interp
view(-23,27)
lighting none

% Here you can see the gaps in the macroalgal surface because we haven't
% patched it in the same way we did the coral one.

%% Collecting bifurcation data
% In the main text figure, we include side panels that show the bifurcation
% diagrams along a specific transect. Here, we generate those model outputs
% for plotting.
refset2 = linspace(0,1,100);

hystchoice = round(length(herbset)/10*3,0);
h = herbset(hystchoice); % update herbivory value

CLstar_high = NaN*refset2;
CLstar_low = CLstar_high;
CLstar_un = CLstar_high;

Mstar_high = NaN*refset2;
Mstar_low = NaN*refset2;
Mstar_un = NaN*refset2;

eqctr = 0;
    for j = 1:length(refset2)
        p = refset2(j); % update refuge value
        
        eq1 = (gC * CL + gammaC) * (1 - CL - CP - CD - ME - MP) - gMC * CL * (ME + MP) - dC * CL==0;
        eq2 = gC*CP*CD + s*gammaC*CD - gMC*CP*(ME+MP) - dC*CP - er*CP == 0;
        eq3 = dC * (CL + CP) + MP * (h * (1 - p) + dM) -  s*gammaC*CD - gC*CP*CD - gM * CD * (ME + MP) - er * CD -  gammaM * CD==0;
        eq4 = (gM * ME + gM * MP + gammaM) * (1 - CL - CP - CD - ME - MP) - ME * (h + dM)==0;
        eq5 = gM * (ME + MP) * CD + gMC * (CL+CP) * (ME + MP) + gammaM * CD - MP * (h * (1 - p) + dM + er)==0;

        S = vpasolve([eq1, eq2, eq3, eq4, eq5],[CL, CP, CD, ME, MP]); % VPASOLVE finds equilibria numerically
        res = [S.CL,S.CP,S.CD,S.ME,S.MP]'; % stores the equilibria in a matrix. We use the transpose here so that the rows are state variables, and the columns are equilibria
        res = res(:,abs(sum(imag(res)))==0); % Remove all imaginary roots
        res = res(:,min(res)>=0); % Remove all negative roots

if size(res,2) == 1
    if eqctr == 1
        CLstar_low(j) = res(1,1)+res(2,1);
        Mstar_low(j) = res(4,1)+res(5,1);
    else
        CLstar_high(j) = res(1,1)+res(2,1);
        Mstar_high(j) = res(4,1)+ res(5,1);
    end
else
    eqctr = 1;
    CLstars = sort(res(1,:)+res(2,:));
    Mstars = res(4,:)+res(5,:);
    Mstars = sort(Mstars);
    CLstar_low(j) = CLstars(1);
    CLstar_un(j) = CLstars(2);
    CLstar_high(j) = CLstars(3);
    Mstar_low(j) = Mstars(3);
    Mstar_un(j) = Mstars(2);
    Mstar_high(j) = Mstars(1);
    
end
        
    end
    
%% Make the full plot

% This set of commands generates a figure of the correct height and width
figure(2)
close(2)
figure(2)
coords = get(gcf,'Position');
set(gcf,'Position',[coords(1),coords(2),coords(3)*1.5,coords(4)])

% Set up plotting detals
herblim = 0.4; % plot herbivory to h = 0.4
corcol = [107/255,126/255,158/255]; % Coral color is blue
maccol = [144/255,126/255,81/255]; % Macroalgal color is yellow

% Plot the coral bifurcation diagram
subplot(2,6,[1:2])
plot(refset2,CLstar_high,'LineWidth',2,'Color',corcol)
hold on
plot(refset2,CLstar_low,'LineWidth',2,'Color',corcol)
plot(refset2,CLstar_un,'LineStyle','--','LineWidth',2,'Color',corcol)
xlabel('Algal Protection (\rho)')
ylabel('Coral (C_L+C_D)')
ylim([0,1])

% Plot the macroalgal bifurcation diagram
subplot(2,6,[7:8])
plot(refset2,Mstar_high,'LineWidth',2,'Color',maccol)
hold on
plot(refset2,Mstar_low,'LineWidth',2,'Color',maccol)
plot(refset2,Mstar_un,'LineStyle','--','LineWidth',2,'Color',maccol)
xlabel('Algal Protection (\rho)')
ylabel('Macroalgae (M_S+M_C)')
ylim([0,1])

% Plot the surface
subplot(2,7,[4:7,11:14])
surf(refset,herbset,highstab_CL','EdgeColor','n','FaceLighting','flat')
hold on
surf(refset,herbset,lowstab_CL')
surf(refset,herbset,unstab_mod_CL','FaceColor','red','EdgeColor','n')
ylabel('Herbivores (h)')
xlabel('Algal Protection (\rho)')
zlabel('Coral Cover (C_L + C_D)')
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
lowhyst = lowstab_CL(:,hystchoice);
lowhyst = lowhyst(isfinite(lowhyst));
lowhystx = refset(isfinite(lowstab_CL(:,hystchoice)));

% select unstable edge
unhyst = unstab_CL(:,hystchoice);
unhyst = unhyst(isfinite(unhyst));
unhystx = refset(isfinite(unstab_CL(:,hystchoice)));

% select upper edge
highhyst = highstab_CL(:,hystchoice);
highhyst = highhyst(isfinite(highhyst));
highhystx = refset(isfinite(highstab_CL(:,hystchoice)));

% concatenate and plot
xset = [highhystx,wrev(unhystx),lowhystx];
yset = herbset(hystchoice)*ones(length(xset));
zset = [highhyst',wrev(unhyst'),lowhyst'];
%plot3(xset,yset,zset,'LineWidth',3,'Color',[0.6350 0.0780 0.1840])
plot3(xset,yset,zset,'LineWidth',3,'Color','k')


% We did some further aesthetic modifications post-plotting in PowerPoint
% to help show the hysteresis line above the surface, and to patch some
% continued holes in the surface.




%% Bifurcation in s
% In the main text figure, we include side panels that show the bifurcation
% diagrams along a specific transect. Here, we generate those model outputs
% for plotting.
gC = 0.1;
dC = 0.001;
gammaC = 0.001;
gM = 0.6;
gMC = 0.06;
dM = 0.05;
h = 0.3;
gammaM = 0.001;
er = 0.01;
p = 0.5;


s_set = linspace(0,5,100);


CLstar_high = NaN*s_set;
CLstar_low = CLstar_high;
CLstar_un = CLstar_high;

Mstar_high = NaN*s_set;
Mstar_low = NaN*s_set;
Mstar_un = NaN*s_set;

eqctr = 0;
    for j = 1:length(s_set)
        s = s_set(j); % update refuge value
        
        eq1 = (gC * CL + gammaC) * (1 - CL - CP - CD - ME - MP) - gMC * CL * (ME + MP) - dC * CL==0;
        eq2 = gC*CP*CD + s*gammaC*CD - gMC*CP*(ME+MP) - dC*CP - er*CP == 0;
        eq3 = dC * (CL + CP) + MP * (h * (1 - p) + dM) -  s*gammaC*CD - gC*CP*CD - gM * CD * (ME + MP) - er * CD -  gammaM * CD==0;
        eq4 = (gM * ME + gM * MP + gammaM) * (1 - CL - CP - CD - ME - MP) - ME * (h + dM)==0;
        eq5 = gM * (ME + MP) * CD + gMC * (CL+CP) * (ME + MP) + gammaM * CD - MP * (h * (1 - p) + dM + er)==0;

        S = vpasolve([eq1, eq2, eq3, eq4, eq5],[CL, CP, CD, ME, MP]); % VPASOLVE finds equilibria numerically
        res = [S.CL,S.CP,S.CD,S.ME,S.MP]'; % stores the equilibria in a matrix. We use the transpose here so that the rows are state variables, and the columns are equilibria
        res = res(:,abs(sum(imag(res)))==0); % Remove all imaginary roots
        res = res(:,min(res)>=0); % Remove all negative roots

if size(res,2) == 1
    if eqctr == 1
        CLstar_low(j) = res(1,1)+res(2,1);
        Mstar_low(j) = res(4,1)+res(5,1);
    else
        CLstar_high(j) = res(1,1)+res(2,1);
        Mstar_high(j) = res(4,1)+ res(5,1);
    end
else
    eqctr = 1;
    CLstars = sort(res(1,:)+res(2,:));
    Mstars = res(4,:)+res(5,:);
    Mstars = sort(Mstars);
    CLstar_low(j) = CLstars(1);
    CLstar_un(j) = CLstars(2);
    CLstar_high(j) = CLstars(3);
    Mstar_low(j) = Mstars(3);
    Mstar_un(j) = Mstars(2);
    Mstar_high(j) = Mstars(1);
    
end
        
    end
    
    
%%
    
figure(3)
close(3)
figure(3)
plot(s_set,CLstar_low)
hold on
plot(s_set,CLstar_un,'LineStyle','--')
plot(s_set,CLstar_high)
    