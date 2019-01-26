% Dibuat oleh: Sigit Prakarsa
% Alumni: T. Aeronotika Astronotika ITB
% Kontak: IG: @sgtprakarsa
% 25 Januari 2019
% Code ini dapat disebarluaskan untuk tujuan pendidikan

%***********************************************************************************%
% COPYRIGHT: 
%       engcode-ID | Coding untuk engineer indonesia
%       github  : https://github.com/engcode-ID
%       IG      : https://www.instagram.com/engcode.id/
%
%***********************************************************************************%

close all;
clear all;
clc;

%***********************************************************************************%
%
% ANALISIS TRUSS 2 DIMENSI MENGGUNAKAN METODE STIFFNESS
% 
% Analisis truss (batang) adalah salah satu analisis struktur 
% sederhana untuk menentukan displacement, deformation, ataupun 
% gaya yang bekerja pada setiap batang.
%
% Ada beberapa metode dalam penyelesaian sistem truss ini, seperti
% metode Force (flexibility) dan Displacement (stiffness).
%
% User akan diminta untuk menginput:
%   1. koordinat truss (nodal dan hubungan antar nodal)
%   2. gaya yang bekerja pada tiap nodal
%   3. kondisi batas (support, dsj)
%   4. properti stiffness (modulus elastis dan cross-section area)
%
% Selanjutnya, user akan dapat melihat:
%   1. Stiffness Matrix
%   2. Displacement vector
%   3. Nodal force vector
%   4. Element force
%   5. Element stress
%   6. Tampilan sistem truss
%
%***********************************************************************************%
%
% REFERENSI
%
% https://engineering.purdue.edu/~aprakas/CE474/CE474-Ch5-StiffnessMethod.pdf
%
%***********************************************************************************%

%***********************************************************************************%
% Input 
%***********************************************************************************%
% Masukkan lokasi joint/nodal dari truss
% Format sebagai berikut
% Node = [Ai Bi Ci;
%         Ai+1 Bi+1 Ci+1;
%         ...]
% Dimana: 
%   Ai = Nodal ke-i, 
%   Bi = koordinat x untuk nodal ke-i, 
%   Ci = koordinat y untuk nodal ke-i
%   Semua satuan dalam m
Node = [1 0 4;
        2 3 4;
        3 0 0;
        4 3 0;
        5 6 0];

% Masukkan konektivitas dari tiap nodal untuk membuat elemen truss
% Satu elemen truss (1 batang) tersusun dari 2 nodal yang dihubungkan.
% Format sebagai berikut 
% Element = [Ai Bi Ci;
%           Ai+1 Bi+1 Ci+1;
%           ...]
% Dimana: 
%   Ai = Nomor elemen truss, 
%   Bi = Nodal pertama untuk elemen ke-i, 
%   Ci = Nodal kedua untuk elemen ke-i
Element =[1 1 2;
          2 1 3;
          3 3 2;
          4 2 4;
          5 2 5;
          6 3 4;
          7 4 5];

% Masukkan gaya yang bekerja pada tiap nodal 
% gaya dimasukkan dengan format sebagai berikut.
% FInput = [Ai Bi;
%           Ai+1 Bi+1;
%           ...]
% dimana: 
%   Ai adalah penunjuk lokasi gaya.
%   Ai = 2*indeks nodal - 1     untuk gaya horizontal
%   Ai = 2*indeks nodal         untuk gaya vertikal
%   Bi = besar gaya yang bekerja (Newton)
%
%   Contoh:
%     apabila ada gaya sebesar 100 N yang bekerja pada nodal ke-3 searah sumbu x, 
%     maka anda masukkan:
%        indeks nodal = 3, gaya = horizontal (searah sumbu-x)
%         Ai = 2*3 - 1 = 5   
%         Bi = 100
%     FInput = [5 100];   
FInput = [3 10000;
          8 -10000];

% Input kondisi batas (displacement)
% Masukkan kondisi batas (fixed support, roll support, dsj) dengan format sebagai berikut
% BCInput = [Ai Bi;
%           Ai+1 Bi+1;
%           ...]
% dimana: 
%   Ai adalah penunjuk lokasi kondisi batas.
%   Ai = 2*indeks nodal - 1     untuk kondisi batas horizontal
%   Ai = 2*indeks nodal         untuk kondisi batas vertikal
%   Bi = 0
%
%   Contoh:
%     apabila ada fixed support pada nodal ke-3, 
%     maka anda masukkan:
%        indeks nodal = 3, 
%        karena fixed support:
%           displacement arah horizontal = 0
%           displacement arah vertikal = 0
%         
%         > arah horizontal
%         Ai = 2*3 - 1 = 5   
%         Bi = 0
%         > arah vertikal
%         Ai = 2*3 = 6   
%         Bi = 0
%     BCInput = [5 0;
%                 6 0];   
BCInput = [1 0;
          2 0;
          5 0;
          9 0;
          10 0];

% Masukan properti stiffness elemen
% properti stiffness didefinisikan sebagai modulus elastis (E) dan cross section area (A) pada
% tiap elemen. Masukkan nilai-nilai tersebut dengan format sebagai berikut
% Prop = [Ai Bi Ci;
%         Ai+1 Bi+1 Ci+1;
%         ...];
% dimana : 
%   Ai = elemen ke-i, 
%   Bi = modulus elasti [Pa] pada elemen ke-i, 
%   Ci = Cross section area [m^2] pada elemen ke-i
Prop =[1 2e11 1e-5;
      2 2e11 1e-5;
      3 2e11 1e-5;
      4 2e11 1e-5;
      5 2e11 1e-5;
      6 2e11 1e-5;
      7 2e11 1e-5];

%***********************************************************************************%
% METODE PENYELESAIAN
% 1) definisikan jumlah DOF (degree of freedom) pada tiap elemen
TotalNode = size(Node, 1); 
TotalElement = size(Element, 1);
ElDOF = zeros(TotalElement, 4);
for i = 1:TotalElement
    ElDOF(i, 1) = 2 * Element(i, 2) - 1;   %DOF horizontal untuk nodal pertama pada elemen ke-i
    ElDOF(i, 2) = 2 * Element(i, 2);       %DOF vertikal untuk nodal pertama pada elemen ke-i

    ElDOF(i, 3) = 2 * Element(i, 3) - 1;   %DOF horizontal untuk nodal kedua pada elemen ke-i
    ElDOF(i, 4) = 2 * Element(i, 3);       %DOF vertikal untuk nodal kedua pada elemen ke-i
end

% 2) definisikan matriks stiffness, gaya, dan perpindahan (displacement)
Force = zeros(2 * TotalNode, 1);
for i = 1:size(FInput, 1)
    Force(FInput(i, 1), 1) = FInput(i, 2);
end 

Displacement = zeros(2 * TotalNode, 1);
for i = 1:size(BCInput, 1)
    Displacement(BCInput(i, 1), 1) = BCInput(i, 2);
end

Stiffness = zeros(2 * TotalNode);

A = Prop(:, 3);
E = Prop(:, 2);

% 3) definisikan koordinat elemen
TrussX = zeros(TotalElement, 2);
TrussY = zeros(TotalElement, 2);

% 4) mencari nilai dari matriks stiffness elemen pada koordinat global
for i = 1:TotalElement
    TrussX(i, 1) = Node(Element(i, 2), 2);
    TrussX(i, 2) = Node(Element(i, 3), 2);
    TrussY(i, 1) = Node(Element(i, 2), 3);
    TrussY(i, 2) = Node(Element(i, 3), 3);

    L(i, 1) = sqrt((TrussX(i, 1) - TrussX(i, 2))^2 + (TrussY(i, 1) - TrussY(i, 2))^2);    % Mencari panjang elemen
    klocal(i, 1) = A(i, 1) * E(i, 1) / L(i, 1);                                           % Mencari nilai matriks stiffness dari elemen pada koordinat lokal
    C(i, 1) = (TrussX(i, 2) - TrussX(i, 1)) / L(i, 1);                                    % Mencari nilai cosine untuk matriks transformasi
    S(i, 1) = (TrussY(i, 2) - TrussY(i, 1)) / L(i, 1);                                    % Mencari nilai sinus untuk matriks transformasi
    
    % membuat matriks transformasi
    T = [C(i,1)^2 C(i,1)*S(i,1) -(C(i,1)^2) -C(i,1)*S(i,1);
        C(i,1)*S(i,1) S(i,1)^2 -C(i,1)*S(i,1) -(S(i,1)^2);
        -(C(i,1)^2) -C(i,1)*S(i,1) C(i,1)^2 C(i,1)*S(i,1);
        -C(i,1)*S(i,1) -(S(i,1)^2) C(i,1)*S(i,1) S(i,1)^2];

    kglobal = klocal(i, 1) * T;                                                              % transformasi matriks dari koodrinat lokal ke global
    for j = 1:4
       for k = 1:4
          Stiffness(ElDOF(i, j), ElDOF(i, k)) = Stiffness(ElDOF(i, j), ElDOF(i, k)) + kglobal(j, k);
       end
    end
end

% 5) menentukan perpindahan pada DOF lain
for i = 1:2*TotalNode
    UnknownDisp(i, 1) = i;
end
UnknownDisp(BCInput(:, 1)) = [];

for i = 1:size(UnknownDisp, 1) 
    disp_new(i, 1) = Displacement(UnknownDisp(i, 1), 1);    % mendefinisikan the undefined DOF Vector
    force_new(i, 1) = Force(UnknownDisp(i, 1), 1);          % mendefinisikan gaya pada DOF tersebut
end 

for i = 1:size(UnknownDisp, 1)
    for j = 1:size(UnknownDisp, 1)
        stiff_new(i, j) = Stiffness(UnknownDisp(i, 1), UnknownDisp(j, 1));   
    end
end
disp_new = stiff_new \ force_new;

% 6) menentukan perpindahan nodal
for i = 1:size(UnknownDisp, 1)
    Displacement(UnknownDisp(i, 1), 1) = disp_new(i, 1);
end  

% 7) menentukan reaksi pada Support
KnownDisp = BCInput(:, 1);
for i = 1:size(KnownDisp, 1)
    Force(KnownDisp(i, 1), 1) = Stiffness(KnownDisp(i, 1), :) * Displacement;
end
 
% 8) Menentukan Stress pada tiap elemen
for i = 1:TotalElement
    varA = [C(i, 1) S(i, 1) -C(i, 1) -S(i, 1)];
    varB = [Displacement(ElDOF(i, 3), :);
            Displacement(ElDOF(i, 4), :);
            Displacement(ElDOF(i, 1), :);
            Displacement(ElDOF(i, 2), :)];
    varC = L(i, 1) * E(i, 1);
    ElStress(i, 1) = varA * varB / varC;
    ElForce(i, 1)= ElStress(i, 1) * A(i, 1);
end

%***********************************************************************************%
% Menampilkan hasil
disp 'Stiffness Matrix'
disp (Stiffness)
disp 'Displacement Vector'
disp (Displacement)
disp 'Nodal Force Vector'
disp (Force)
disp 'Element Force'
disp (ElForce)
disp 'Element Stress'
disp (ElStress) 

% Gambar
figure(1);
hold on;
plot(TrussX, TrussY, 'k.')
hold on;
axis equal;
for el = 1:TotalElement
    elnodes = Element(el, 2:3);
    nodexy = Node(elnodes, :);
    plot(nodexy(:, 2), nodexy(:, 3), 'k--')
end;

% gambar baru
%magni = 20;
%nodeNew = node + magni*reshape(Displacement, 2, TotalNode)
%disp(nodeNew)
magnification = 20;
newDisplacement = reshape(Displacement, 2, TotalNode)
for i = 1:TotalNode
    nodeNew(i, 1) = Node(i, 1);
    nodeNew(i, 2) = Node(i, 2) + magnification*newDisplacement(1, i);
    nodeNew(i, 3) = Node(i, 3) + magnification*newDisplacement(2, i);
end;

plot(nodeNew(:, 2), nodeNew(:, 3), 'o')
hold on;
axis equal;
for el = 1:TotalElement
    elnodes = Element(el, 2:3);
    nodexy = nodeNew(elnodes, :);
    plot(nodexy(:, 2), nodexy(:, 3), 'k-')
end
