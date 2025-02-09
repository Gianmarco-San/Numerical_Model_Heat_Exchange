function [D,bconst, bvar] = BuildStiffEvol_cil(Me)
%Assemble the matrix D and the vector b of the Diffusion problem with
%homogeneous B.C.s
%Input:
%   Me     :a Mesh2D object
%
%Output:
%   D      :diffusion matrix
%   b      :constant terms vector

%for clarity, call some properties of Me with shorter names
V=Me.Triangles.Vertices;
Areas=Me.Triangles.Areas;
CenterOfMass=Me.Triangles.CenterOfMass;
Nodes=Me.Nodes;
Dof=Me.Nodes.Dof;
%number of internal nodes: we know that the N unknown nodes are numbered from
%1 to N in Me.UnknownNodes; the maximum is therefore the number of unknown
%(degrees of freedom)
numDof = max(Dof);

%vectors preallocation: instead of allocating the (sparse) diffusion matrix, 
%we save the rows, columns and values corresponding to each contribution; 
%at the end, we'll call sparse(...) to obtain the diffusion matrix
bconst = zeros(numDof,1);
bvar = zeros(numDof,1);
row = zeros(Me.MatrixContributions,1);
col = zeros(Me.MatrixContributions,1);
d = zeros(Me.MatrixContributions,1);
pos=1;  %we start from the element in position 1, we'll increase this index 
        %everytime we add an entry

mu=Me.mu;
beta=Me.beta;
rho=Me.rho;
%main loop on each triangle 

for e=1:size(V,1)   
    Dz(1) = Nodes.Y(V(e,3)) - Nodes.Y(V(e,2));
    Dz(2) = Nodes.Y(V(e,1)) - Nodes.Y(V(e,3));
    Dz(3) = Nodes.Y(V(e,2)) - Nodes.Y(V(e,1));
    Dr(1) = Nodes.X(V(e,3)) - Nodes.X(V(e,2));
    Dr(2) = Nodes.X(V(e,1)) - Nodes.X(V(e,3));
    Dr(3) = Nodes.X(V(e,2)) - Nodes.X(V(e,1));
    
    %we evaluate the external force in the center of mass of this triangle
    
    rb=CenterOfMass.X(e);   

    %for each vertex of this triangle 
    for ni=1:3
        %look at the "unknown" numbering: if the node is positive, it
        %corresponds to a degree of freedom of the problem
        ii = Dof(V(e,ni));
        %is it unknown?
        if ii > 0 
            %yes it is! second loop on the vertices
            for nj=1:3    
                jj = Dof(V(e,nj));               
                betar=beta(e,1)-mu(e)/rb;
                betay=beta(e,2);
                diffusion=mu(e)*(Dr(ni)*Dr(nj)+Dz(ni)*Dz(nj))/(4.0*Areas(e));
                transport=(betar*Dz(nj)-betay*Dr(nj))*rho(e)/6;
                 %%is it unknown as well?
                if jj > 0                                   
                    %Non sparse solution: D(ii,jj)=D(ii,jj) +  diffusion + transport;
                    row(pos)=ii;
                    col(pos)=jj;                    
                    d(pos)= diffusion + transport;
                    pos=pos+1;
                else
                    value=Me.BC.DirichletNodes(-jj,2);
                        bvar(ii) = bvar(ii) - (diffusion+transport)*value;
                end
            end
            
        end
    end   
end



% Devo aggiungere questo ciclo for per tener conto delle condizioni di
% Robin
Edges=Me.Edges;         % Matrice di due colonne con i vertici dei lati della mesh
Robin=Me.BC.RobinEdges; % Matrice lunga come i lati di Robin: la prima colonna tiene l'indice di lato, e le altre due tengono i parametri della condizione di Robin nell'ordine assegnato

for k=1:size(Robin,1)           %dopo aver fatto il ciclo sui triangoli faccio il ciclo sui lati di robin
    Node1=Edges(Robin(k,1),1);  %Inizio Calcolo della lunghezza del lato: estraggo i nodi,  
    Node2=Edges(Robin(k,1),2);
    dx=Nodes.X(Node1)-Nodes.X(Node2);
    dy=Nodes.Y(Node1)-Nodes.Y(Node2);    
    dist=sqrt(dx*dx+dy*dy);     % Fine calcolo della lunghezza del lato
    ii1=Dof(Node1); 
    ii2=Dof(Node2);
    g=Robin(k,3);
    h=Robin(k,2); % Vedo dove sta il lato di Robin: devo vedere se tutti e due gli estremi sono gradi di libertà o no:
    if ii1>0 && ii2<0 %ii1 is unknown, ii2 is known
%         Dirichlet=Me.BC.Dirichlet(-ii2,2); %Devo tenere conto dell'interazione tra condizioni di Robin e Dirichlet non omogeneo (non nel nostro caso)
        bconst(ii1)=bconst(ii1)+g/2*dist-h*dist/6*Dirichlet;
        row(pos)=ii1;
        col(pos)=ii1;
        d(pos)=h*dist/3;
        pos=pos+1;
        %D(ii1,ii1)=D(ii1,ii1)+h*dist/3;
    elseif ii1<0 && ii2>0 %ii1 is known, ii2 is unknown
%         Dirichlet=Me.BC.Dirichlet(-ii1,2); %Devo tenere conto dell'interazione tra condizioni di Robin e Dirichlet non omogeneo (non nel nostro caso)
        bconst(ii2)=bconst(ii2)+g/2*dist-h*dist/6*Dirichlet;        
        row(pos)=ii2;
        col(pos)=ii2;
        d(pos)=h*dist/3;
        pos=pos+1;
        %D(ii2,ii2)=D(ii2,ii2)+h*dist/3;
    else  %both are unknwon
        bconst(ii1)=bconst(ii1)+g/2*dist;
        bconst(ii2)=bconst(ii2)+g/2*dist;
        row(pos:pos+3)=[ii1;ii2;ii1;ii2];
        col(pos:pos+3)=[ii1;ii2;ii2;ii1];
        d(pos:pos+3)=[2;2;1;1]*h*dist/6;
        pos=pos+4;
        %D(ii1,ii1)=D(ii1,ii1)+h*dist/3;
        %D(ii2,ii2)=D(ii2,ii2)+h*dist/3;
        %D(ii1,ii2)=D(ii1,ii2)+h*dist/6;
        %D(ii2,ii1)=D(ii2,ii1)+h*dist/6;        
    end
end
%assemble the stiffness matrix D from the

D=sparse(row,col, d, numDof, numDof);