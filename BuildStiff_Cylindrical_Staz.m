function [D,b] = BuildStiff_Cylindrical_Staz(Me)
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
numDof = max(Dof); % Numero di gradi di libertà

%vectors preallocation: instead of allocating the (sparse) diffusion matrix, 
%we save the rows, columns and values corresponding to each contribution; 
%at the end, we'll call sparse(...) to obtain the diffusion matrix
b = zeros(numDof,1);
row = zeros(Me.MatrixContributions,1);
col = zeros(Me.MatrixContributions,1);
d = zeros(Me.MatrixContributions,1);
pos=1;  %we start from the element in position 1, we'll increase this index 
        %everytime we add an entry

mu =Me.mu;
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
                betar = beta(e,1)-mu(e)/rb; 
                betaz = beta(e,2);

                diffusion = mu(e)*(Dr(ni)*Dr(nj)+Dz(ni)*Dz(nj))/(4.0*Areas(e));
                transport = - (betaz*Dr(nj)-betar*Dz(nj))*rho(e)/6; 
                
                 %%is it unknown as well?
                if jj > 0                                   
                    row(pos)=ii;
                    col(pos)=jj;                    
                    d(pos)= diffusion + transport;
                    pos=pos+1;
                else
                    value=Me.BC.DirichletNodes(-jj,2); 
                    b(ii)=b(ii)-(diffusion+transport)*value;
                end
            end
        end
    end
end

%% Robin B.C.s
Edges=Me.Edges; % Mi dice i vertici di ogni segmento
Robin=Me.BC.RobinEdges; % E' una matrice di 3 colonne (1. segmento; 2. coeff. conv.; 3. h*T)
for k=1:size(Robin,1)
    Node1=Edges(Robin(k,1),1);
    Node2=Edges(Robin(k,1),2);
    dx=Nodes.X(Node1)-Nodes.X(Node2);
    dy=Nodes.Y(Node1)-Nodes.Y(Node2);    
    dist=sqrt(dx*dx+dy*dy);
    ii1=Dof(Node1);
    ii2=Dof(Node2);
    g=Robin(k,3);
    h=Robin(k,2);
    if ii1>0 && ii2<0 %ii1 is unknown, ii2 is known
        
%        Dirichlet = Me.BC.Dirichlet(-ii2,2); %aggiunto per quando ci sono nodi di Dirichlet insieme a Robin
        b(ii1)=b(ii1)+g/2*dist-h*dist/6*Dirichlet;
        row(pos)=ii1;
        col(pos)=ii1;
        d(pos)=h*dist/3;
        pos=pos+1;
        %D(ii1,ii1)=D(ii1,ii1)+h*dist/3;
    elseif ii1<0 && ii2>0 %ii1 is known, ii2 is unknown
        
%        Dirichlet = Me.BC.Dirichlet(-ii1,2); %aggiunto per quando ci sono nodi di Dirichlet insieme a Robin
        b(ii2)=b(ii2)+g/2*dist-h*dist/6*Dirichlet;      
        
        row(pos)=ii2;
        col(pos)=ii2;
        d(pos)=h*dist/3;
        pos=pos+1;
  
    else  %both are unknwon
        b(ii1)=b(ii1)+g/2*dist;
        b(ii2)=b(ii2)+g/2*dist;
        row(pos:pos+3)=[ii1;ii2;ii1;ii2];
        col(pos:pos+3)=[ii1;ii2;ii2;ii1];
        d(pos:pos+3)=[2;2;1;1]*h*dist/6;
        pos=pos+4;
            
    end
end

%assemble the stiffness matrix D from the
D=sparse(row,col, d, numDof, numDof);