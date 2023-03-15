function Programme_Principal_v2

%%% Affichage du petit cadre où l'on entre nos données%%%
prompt={'Entrer le nom du fichier accélération:','Entrer le nomdu fichier gyroscope:','Entrer la fréquence d échantillonage:'};
dlgtitle='données d entrée';
dims=[1 35];
definput={'S8_T2_14_Thorax_acc.csv','S8_T2_14_Gyro.csv','208'};%%Donne les valeurs par défaut, ici les noms des fichiers tests
answer=inputdlg(prompt,dlgtitle,dims,definput);%liste recevant toutes les données entrées par l'utilisateur
nom_fichier_acc=cell2mat(answer(1));%%passage de la variable 1 de Answer (pour le nom du fichier accélération) en une matrice
nom_fichier_gyr=cell2mat(answer(2));%%passage de la variable 2 de Answer (pour le nom du fichier gyroscope) en une matrice
freq=str2double(cell2mat(answer(3)));%%passage de la variable 3 de Answer (pour la fréquence) en une matrice puis passage d'une chaine de carractères en num (double est plus performant)

%%ouverture du fichier demandé en format csv
Acc=readtable(nom_fichier_acc); %%ouverture du fichier accélération
Gyro=readtable(nom_fichier_gyr); %%ouverture du fichier gyroscope

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%définition des variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Variables accélérations
AY=Acc.y;%%on récupère seulement les valeurs de l'accélération suivant y
AZ=Acc.z;%%on récupère seulement les valeurs de l'accélération suivant z

%%mise en place de la base de temps
t=zeros([1 length(AY)]);
j=1;
for i=2:length(t)
    t(i)=j*(1/freq);
    j=j+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%filtration des signaux%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fc=10;
[b,a]=butter(4,fc/(freq/2));
AyF=filtfilt(b,a,AY); %%signal vertical filtré avec zero-lag

[b,a]=butter(4,fc/(freq/2));
AzF=filtfilt(b,a,AZ); %%signal anterior-posterior filtré avec zero-lag

%%on ajuste la base de temps et les courbes pour que ça colle avec le début de
%%l'expérience%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start=4162;
taj=t(start:length(AyF));
Ayaj=AY(start:end);
AyFaj=AyF(start:end);
AzFaj=AzF(start:end);
Azaj=AZ(start:end);
G=Gyro.z(start:length(Ayaj)+start-1);%%Cela permet à G et A d'avoir la même taille

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% recherches des points carractéristiques %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%On cherche les pics d'accélération verticale
[~,Maxy]=findpeaks(AyFaj,'MinPeakHeight',10);

%%On inverse la courbe
AzFajInv=inverser_fonction(AzFaj);
AyFajInv=inverser_fonction(AyFaj);

%%On cherche les minimums locaux (qui sont les maximums de la fonction
%%inverse)
[~,Minz,~,~]=findpeaks(AzFajInv);
[~,Miny]=findpeaks(AyFajInv,'MinPeakHeight',1);

%%% Tri dans les points de décollage et d'amortissement entre deux pics de
%%% Maxy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Error=[];
decollage=[];
amortissement=[];
NMaxy=[];
Maxz=[];
DecollageSolo=[];
for k=1:length(Maxy)-2
    index=[];
    index2=[];
    diff=[];
    Num=[];
    for l=1:length(Minz)-1
        if (Minz(l)>=Maxy(k))&&(Minz(l)<=Maxy(k+1))
            index(end+1)=Minz(l);
            Num(end+1)=l;
        end
    end
    for m=1:length(Miny)
        if (Miny(m)>=Maxy(k))&&(Miny(m)<Maxy(k+1))
            index2(end+1)=Miny(m);
        end
    end
    if (length(index)>=2)&&(length(index2)>=1)
        NMaxy(end+1)=index2(1);
        for n=1:length(index)
            diff(end+1)=abs(index(n)-index2(1));
        end
        [~,p]=min(diff);
        decollage(end+1)=index(p);
        if p<length(index)
            amortissement(end+1)=index(p+1);
        else
            Error(end+1)=Maxy(k);
            Error(end+1)=Maxy(k+1);
            DecollageSolo(end+1)=index(p);
            NA=AzFaj(index(p):Minz(Num(length(Num))+1));
            [~,IndexMax]=max(NA);
            Maxz(end+1)=IndexMax(1)+index(p);
        end
    elseif (length(index)==1)&&(length(index2)>=1)
        decollage(end+1)=index(1);
        DecollageSolo(end+1)=index(1);
        Error(end+1)=Maxy(k);
        Error(end+1)=Maxy(k+1);
        NA=AzFaj(index(1):Minz(Num(length(Num))+1));
        [~,IndexMax]=max(NA);
        Maxz(end+1)=IndexMax(1)+index(1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Méthode altérnative pour trouver le minimum %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AzajInv=inverser_fonction(Azaj);
[~,MinzNf,~,~]=findpeaks(AzajInv);


PositionMini=zeros(1,length(DecollageSolo));
for x=1:length(DecollageSolo)
    taille=Maxz(x)-DecollageSolo(x)-8;
    Interv=zeros(1,taille);
    M=1;
    c=1;
    for y=DecollageSolo(x)+3:Maxz(x)-4
        z=AzFaj(y+1)-AzFaj(y);
        Interv(M)=z;
        M=M+1;
    end
    while Interv(c)<Interv(c+1)
        c=c+1;
        if c==taille-2
            break;
        end
    end
    while Interv(c)>Interv(c+1)
        c=c+1;
        PositionMini(x)=c+DecollageSolo(x)+3;
        if c==taille-1
            break;
        end
    end
end

%%%Affinage des points en corrélations avec les minimums locaux déterminés
%%%précédemment%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PositionMiniAj=zeros(1,length(PositionMini));
for e=1:length(PositionMini)
    Interv2=[];
    for e2=1:length(MinzNf)
        if (MinzNf(e2)>=PositionMini(e))&&(MinzNf(e2)<Maxz(e))
            Interv2(end+1)=MinzNf(e2);
        end
    end
    if length(Interv2)>=1
        PositionMiniAj(e)=Interv2(1);
    else
        PositionMiniAj(e)=PositionMini(e);
    end
end

%%On le rajoute dans la liste des amortissements %%%%%%%%%%%%%%%%%%%%%%%%%
for f=1:length(PositionMiniAj)
    amortissement(end+1)=PositionMiniAj(f);
end
amortissement=sort(amortissement);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Calcul des temps de vol et de contact %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if decollage(1)<amortissement(1)
    TempsVol=zeros(1,length(amortissement));
    TempsContact=zeros(1,length(amortissement)-1);
    for k5=1:length(TempsVol)-1
        TempsContact(k5)=taj(decollage(k5+1))-taj(amortissement(k5));
    end
    for k6=1:length(TempsVol)
         TempsVol(k6)=taj(amortissement(k6))-taj(decollage(k6));
    end
else
    TempsContact=zeros(1,length(amortissement));
    TempsVol=zeros(1,length(amortissement)-1);
    for k7=1:length(TempsContact)-1
        TempsVol(k7)=taj(amortissement(k7+1))-taj(decollage(k7));
    end
    for k8=1:length(TempsContact)
        TempsContact(k8)=taj(decollage(k8))-taj(amortissement(k8));
    end
end


%%%%%%%%%S'avoir s'il s'agit du pied droit ou du pied gauche %%%%%%%%%%%%%%

Temps_de_Vol_Droit=[];
Temps_de_Contact_Droit=[];
Temps_de_Vol_Gauche=[];
Temps_de_Contact_Gauche=[];
if decollage(1)<amortissement(1)
     xInit=decollage(1);
     xFinal=decollage(2);
      M=moyenne(G(xInit:xFinal));
      if M>=0
          r=1;
      else
          r=2;
      end
      for k9=1:length(TempsVol)-1
          r=r+1;
          if mod(r,2)==0
              Temps_de_Vol_Droit(end+1)=TempsVol(k9);
              Temps_de_Contact_Droit(end+1)=TempsContact(k9);
          else
              Temps_de_Vol_Gauche(end+1)=TempsVol(k9);
              Temps_de_Contact_Gauche(end+1)=TempsContact(k9);
          end
      end
else
    xInit=amortissement(1);
    xFinal=amortissement(2);
    M=moyenne(G(xInit:xFinal));
    if M>=0
       r=1;
    else
       r=2;
    end
    for k10=1:length(TempsContact)-1
          r=r+1;
          if mod(r,2)==0
              Temps_de_Vol_Droit(end+1)=TempsVol(k10);
              Temps_de_Contact_Droit(end+1)=TempsContact(k10);
          else
              Temps_de_Vol_Gauche(end+1)=TempsVol(k10);
              Temps_de_Contact_Gauche(end+1)=TempsContact(k10);
          end
      end
end

    

    function affichage_Vitesse_v2(~,TempsVol,TempsContact,Temps_de_Vol_Droit,Temps_de_Contact_Droit,Temps_de_Vol_Gauche,Temps_de_Contact_Gauche)
        %%on récupère les valeurs calculées grâce à la fonction determination_temps
        %%et on fait les moyennes de toutes ses valeurs
        temps_vol_moyen=moyenne(TempsVol);
        temps_contact_moyen=moyenne(TempsContact);
        TContactPDM=moyenne(Temps_de_Contact_Droit);
        TContactPGM=moyenne(Temps_de_Contact_Gauche);
        TVolPDM=moyenne(Temps_de_Vol_Droit);
        TVolPGM=moyenne(Temps_de_Vol_Gauche);

        %%On calcul le pourcentage de différence entre les valeurs de temps du pied
        %%droit et du pied gauche
        DiffVol=(abs(TVolPDM-TVolPGM)*100)/temps_vol_moyen;
        DiffContact=(abs(TContactPDM-TContactPGM)*100)/temps_contact_moyen;

        %%Création d'un tableau de 3 colonnes et 4 lignes pour contenir toutes les
        %%valeurs calculées.
        temps=zeros(4,2);
        temps(1,1)=temps_vol_moyen;
        temps(1,2)=temps_contact_moyen;

        temps(2,1)=TVolPGM;
        temps(2,2)=TContactPGM;

        temps(3,1)=TVolPDM;
        temps(3,2)=TContactPDM;

        temps(4,1)=DiffVol;
        temps(4,2)=DiffContact;

        %%Affichage de ce tableau sur une nouvelle figure
        fig2=uifigure('Name','Temps caractéristiques','Position',[100, 100, 700, 250]);
        uit=uitable(fig2,'Data',temps);
        uit.ColumnName={'Temps de vol moyen (s)','Temps de contact moyen (s)','Temps de foulée (s)'};
        uit.RowName={'Les deux pieds','Pied gauche','Pied droit','Différence entre les deux pieds (%)'};
        uit.Position=[10, 50, 665, 200];
    end

%%Initialisation des tableaux pour accueillir le nombre de pas gauche et
%%droit
NbPasDroit=zeros(1,1);
NbPasGauche=zeros(1,1);
NbPasDroit(1)=length(Temps_de_Contact_Droit);
NbPasGauche(1)=length(Temps_de_Contact_Gauche);

    function plot_Button_pushed_v2(~,ax,~,taj,decollage,amortissement,G)
        Base_de_Temps_Droit=[];
        Courbe_Droit=[];
        Base_de_Temps_Gauche=[];
        Courbe_Gauche=[];
        d=0;

%%Les lignes qui vont suivre sont là pour déterminer si les valeurs de A
%%correspondent au pied droit ou au pied gauche. On prend le même principe
%%que pour le tri des temps de contact et de temps de vol.

% %%Premier cas de figure, la courbe commence par un temps de contact.

if decollage(1)>amortissement(1)
    xInit2=amortissement(1);
    xFinal2=amortissement(2);
    M=moyenne(G(xInit2:xFinal2));
    if M>=0
        d=1;
    else
        d=0;
    end
    for k13=2:min([length(decollage) length(amortissement)])-1
        xInit2=amortissement(k13);
        xFinal2=amortissement(k13+1);
        d=d+1;
        if mod(d,2)==0
            for m1=xInit2:xFinal2
                Base_de_Temps_Droit(end+1)=taj(m1);
                Courbe_Droit(end+1)=AzFaj(m1);
                Courbe_Gauche(end+1)=NaN; %%Le fait de placer des NaN (Not a Number) permet de laisser un vide dans le tracé
% Si nous n'avions pas mis un NaN, la valeur serait automatiquement 0 
                Base_de_Temps_Gauche(end+1)=NaN;
            end
        else
            for m1=xInit2:xFinal2
                Base_de_Temps_Gauche(end+1)=taj(m1);
                Courbe_Gauche(end+1)=AzFaj(m1);
                Courbe_Droit(end+1)=NaN;
                Base_de_Temps_Droit(end+1)=NaN;
            end
        end
    end
else
    %%%Deuxième cas de figure la courbe commence par un temps de vol.
      xInit2=amortissement(1);
      xFinal2=amortissement(2);
      M=moyenne(G(xInit2:xFinal2));
      if M>=0
          d=1;
      else
          d=0;
      end
    for k14=2:min([length(decollage) length(amortissement)])-1
         d=d+1;
        if mod(d,2)==0
            xInit2=amortissement(k14);
            xFinal2=amortissement(k14+1);
            for m1=xInit2:xFinal2
                Base_de_Temps_Droit(end+1)=taj(m1);
                Courbe_Droit(end+1)=AzFaj(m1);
                Courbe_Gauche(end+1)=NaN;
                Base_de_Temps_Gauche(end+1)=NaN;
            end
        else
            xInit2=amortissement(k14);
            xFinal2=amortissement(k14+1);
             for m1=xInit2:xFinal2
                Base_de_Temps_Gauche(end+1)=taj(m1);
                Courbe_Gauche(end+1)=AzFaj(m1);
                Courbe_Droit(end+1)=NaN;
                Base_de_Temps_Droit(end+1)=NaN;
             end
        end
    end
end

%%On procède à l'affichage sur l'axe créé dans le programme principal
plot(ax,Base_de_Temps_Droit,Courbe_Droit,'g',Base_de_Temps_Gauche,Courbe_Gauche,'r',taj(decollage),AzFaj(decollage),'*r',taj(amortissement),AzFaj(amortissement),'*b')
legend(ax,'Pied droit', 'Pied gauche','Points de poussée', 'Points d amortissement')


    end

%Create a figure window. Pour l'interface utilisateur
fig=uifigure('Name','Décomposition de la course à pied','Unit','centimeters','Position',[5 5 30 30]);

%Creation d'un tableau pour l'affichage du nombre de pas
figPasDroit=uitable('Parent',fig,'Units','centimeters','Position',[23,23,6,2],'Data',NbPasDroit);
figPasDroit.ColumnName={'Nombre de pas droit'};
figPasGauche=uitable('Parent',fig,'Units','centimeters','Position',[23,19,6,2],'Data',NbPasGauche);
figPasGauche.ColumnName={'Nombre de pas gauche'};

%Create a UI axes
ax=uiaxes('Parent',fig,'Units','centimeters','Position',[1,1,22,15]);

%Create a push button
btn=uibutton(fig,'push','Position',[900, 250,100,22],'ButtonPushedFcn',@(btn,event)plot_Button_pushed_v2(btn,ax,AyFaj,taj,decollage,amortissement,G));
btn.Text='Graphique';%%on nomme le bouton (par défaut, c'est button mais c'est pas trop explicite)

fig.CloseRequestFcn=@(src,event)my_closereq(src);%%message qui s'affiche lorsque l'on veut quitter l'interface graphique

%Création d'un nouveau bouton pour avoir les valeurs des vitesses
btn_vitesse=uibutton(fig,'push','Position',[900,190,133,22],'ButtonPushedFcn',@(btn,event)affichage_Vitesse_v2(btn,TempsVol,TempsContact,Temps_de_Vol_Droit,Temps_de_Contact_Droit,Temps_de_Vol_Gauche,Temps_de_Contact_Gauche));
btn_vitesse.Text='Temps caractéristiques';

end









