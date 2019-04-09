function varargout = minips(varargin)
% MINIPS MATLAB code for minips.fig
%      MINIPS, by itself, creates a new MINIPS or raises the existing
%      singleton*.
%
%      H = MINIPS returns the handle to a new MINIPS or the handle to
%      the existing singleton*.
%
%      MINIPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MINIPS.M with the given input arguments.
%
%      MINIPS('Property','Value',...) creates a new MINIPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before minips_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to minips_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help minips

% Last Modified by GUIDE v2.5 08-Apr-2019 20:28:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @minips_OpeningFcn, ...
                   'gui_OutputFcn',  @minips_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before minips is made visible.
function minips_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to minips (see VARARGIN)

% Choose default command line output for minips
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes minips wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = minips_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function namafile_Callback(hObject, eventdata, handles)
% hObject    handle to namafile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of namafile as text
%        str2double(get(hObject,'String')) returns contents of namafile as a double


% --- Executes during object creation, after setting all properties.
function namafile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to namafile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cari.
function cari_Callback(hObject, eventdata, handles)
% hObject    handle to cari (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname] = uigetfile('*.jpg','pilih file');
if isequal([filename,pathname],[0,0])
    return
else 
    fullpath = fullfile(pathname,filename);
    handles.gui.fullpath = fullpath;
    handles.gui.filename = filename;
    handles.gui.pathname = pathname;
    
    set(handles.namafile,'String',handles.gui.fullpath);
    guidata(hObject,handles);
end
% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I;
addpath(handles.gui.pathname);
I = imread(handles.gui.filename);
handles.gui.I = I;
imshow(I)
guidata(hObject,handles);


% --- Executes on slider movement.
function brightness_Callback(hObject, eventdata, handles)
% hObject    handle to brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% get slider value
brightness_value = get(hObject,'Value');

% process
I_edit = handles.gui.I;
if not(islogical(I_edit))
    bright_I(:,:,1) = uint8( bound( double(I_edit(:,:,1)) + brightness_value ) ) ;
    if size(I_edit,3) == 3
        bright_I(:,:,2) = uint8( bound( double(I_edit(:,:,2)) + brightness_value ) ) ;
        bright_I(:,:,3) = uint8( bound( double(I_edit(:,:,3)) + brightness_value ) ) ;
    end
    handles.gui.bright_I = bright_I;
    imshow(bright_I);
    guidata(hObject, handles);
else
    return
end
% --- Executes during object creation, after setting all properties.
function brightness_CreateFcn(hObject, ~, handles)
% hObject    handle to brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in bright_ok.
function bright_ok_Callback(hObject, eventdata, handles)
% hObject    handle to bright_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.brightness,'Value',0);
handles.gui.I = handles.gui.bright_I;
guidata(hObject, handles);

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in simpan.
function simpan_Callback(hObject, eventdata, handles)
% hObject    handle to simpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName] = uiputfile('*.jpg', 'Save As');
if PathName==0
return;
end
Name = fullfile(PathName, FileName);
imwrite(handles.gui.I, Name, 'jpg');
guidata(hObject, handles);

% --- Executes on button press in contrast_ok.
function contrast_ok_Callback(hObject, eventdata, handles)
% hObject    handle to contrast_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.contrast,'Value',0);
handles.gui.I = handles.gui.contrast_I;
guidata(hObject, handles);

% --- Executes on slider movement.
function contrast_Callback(hObject, eventdata, handles)
% hObject    handle to contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% get slider value
contrast_value = get(hObject,'Value'); 
% process 
if contrast_value ~= 0     
    cv = contrast_value / 2;     
    contrast_value = 127 ^ (cv / 100);     
    I_edit = handles.gui.I;    
 if not(islogical(I_edit)) %menyeleksi apakah image bukan biner, eksekusi grayscale    
    contrast_I(:,:,1) = uint8( bound( round( (double(I_edit(:,:,1))-128)*contrast_value + 128 ) ) );     
    if size(I_edit,3)== 3 %jika RGB maka dilanjut
    contrast_I(:,:,2) = uint8( bound( round( (double(I_edit(:,:,2))-128)*contrast_value + 128 ) ) );
    contrast_I(:,:,3) = uint8( bound( round( (double(I_edit(:,:,3))-128)*contrast_value + 128 ) ) ); 
    end 
        handles.gui.contrast_I = contrast_I; 
        imshow(contrast_I); 
        guidata(hObject, handles);
 else
     return
 end
end
% --- Executes during object creation, after setting all properties.
function contrast_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in biner.
function biner_Callback(hObject, eventdata, handles)
% hObject    handle to biner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global threshold; %variable global Thresshold
global gray;%variable global Grayscale
global biner_image; %variable global Biner
edit_biner = handles.gui.I; 
[tinggi,lebar] = size (edit_biner); 
edit_biner = zeros(tinggi,lebar) 
for rows= 1 : tinggi
    for columns= 1 : lebar
        if gray(rows, columns) < threshold; 
           edit_biner (rows, columns) = 0;
        else
    edit_biner (rows, columns) = 1;
        end;
    end;
end;
biner_image=edit_biner;
edit_biner = logical(edit_biner);
imshow(edit_biner);
handles.gui.I = edit_biner;
guidata(hObject, handles);

% --- Executes on button press in grayscale.
function grayscale_Callback(hObject, eventdata, handles)
% hObject    handle to grayscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global threshold; %variable global tresshold
global gray;%variale global grayscale
edite_gray = handles.gui.I; 
try
  [rows, columns, numberOfColorChannels]= size(edite_gray);
  if numberOfColorChannels  == 3
      % convert grayscale
      redChannel = edite_gray(:, :, 1); 
      greenChannel = edite_gray(:, :, 2);
      blueChannel = edite_gray(:, :, 3);
      % Do the weighted sum.
      grayImage = .299*double(redChannel) + ...
                  .587*double(greenChannel) + ...
                  .114*double(blueChannel);
      % You probably want to convert it back to uint8 so you can display it.
      grayImage = uint8(grayImage);
  else
      % It's already gray scale.
      grayImage = edite_gray;  
  end
gray = grayImage;
threshold = graythresh(grayImage)*256;
handles.gui.grayImage = grayImage;
handles.gui.I = handles.gui.grayImage;
imshow(handles.gui.I);
guidata(hObject, handles);
end


% --- Executes on button press in flip_vertikal.
function flip_vertikal_Callback(hObject, eventdata, handles)
% hObject    handle to flip_vertikal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit_flipV = handles.gui.I; 
d=edit_flipV; 
[rows,columns]=size(edit_flipV);  
for k=1:rows
    for b=1:columns
d(rows+1-k,b)=edit_flipV(k,b); 
end
end
handles.gui.I = d;
imshow(d);
guidata(hObject, handles);

% --- Executes on button press in rotate.
function rotate_Callback(hObject, eventdata, handles)
% hObject    handle to rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit_rotate = handles.gui.I; 
[rows,columns,layer]=size(edit_rotate); 
if islogical(edit_rotate)
    d= logical(zeros(columns,rows)) %image biner
else
    d=uint8(zeros(columns,rows,layer)); %image rgb
end
k = rows;
for i=1:rows
    for j=1:columns
        d(j,k,:)=edit_rotate(i,j,:); 
    end
    k=k-1;
end
handles.gui.d = d;
handles.gui.I = handles.gui.d;
imshow(handles.gui.I);
guidata(hObject, handles);


% --- Executes on button press in flip_horizontal.
function flip_horizontal_Callback(hObject, eventdata, handles)
% hObject    handle to flip_horizontal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit_flipH = handles.gui.I;
flip_horizontal =edit_flipH(:,end:-1:1,:);
imshow(flip_horizontal);
handles.gui.I = flip_horizontal;
guidata(hObject, handles);


% --- Executes on button press in dilasi.
function dilasi_Callback(hObject, eventdata, handles)
% hObject    handle to dilasi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

edit_dilasi=handles.gui.I;
kernel= [1 1 1 1 1 1 1;]; %menggunakan titik tengah 1 pada karnel
C=padarray(edit_dilasi,[0 3]); %untuk mengganjal masukan citra secara explisit
Dilasi=false(size(edit_dilasi)); 
for i=1:size(C,1)
    for j=1:size(C,2)-6
        Dilasi(i,j)=sum(kernel&C(i,j:j+6));%menyamakan jumlah karnel tersebut ke dilasi
    end
end
handles.gui.Dilasi = Dilasi;
handles.gui.I = handles.gui.Dilasi;
imshow(handles.gui.I);
guidata(hObject, handles);



% --- Executes on button press in zoom_out.
function zoom_out_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit_zoomout = handles.gui.I;
[rows,columns,layer]=size(edit_zoomout); 
zoomout=uint8(zeros(rows,columns,layer)); 

for i= 1:rows/2
    for j= 1:columns/2
        a = i*2;
        b = j*2;        
        zoomout(i,j,:)=edit_zoomout(a,b,:);
    end
end

handles.gui.zoomout = zoomout;
handles.gui.I = handles.gui.zoomout;
imshow(handles.gui.I);
guidata(hObject, handles);
% --- Executes on button press in zoom_in.
function zoom_in_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit_zoomin = handles.gui.I;
[rows,columns,layer]=size(edit_zoomin); 
zoomin=uint8(zeros(rows,columns,layer)); 
for i= 1:2*rows
    for j= 1:2*columns
        a = i - floor(i/2);
        b = j - floor(j/2);        
        zoomin(i,j,:)=edit_zoomin(a,b,:);
    end
end

handles.gui.zoomin = zoomin;
handles.gui.I = handles.gui.zoomin;
imshow(handles.gui.I);
guidata(hObject, handles);
% --- Executes on button press in median.
function median_Callback(hObject, eventdata, handles)
% hObject    handle to median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pros = handles.gui.I;
[rows,columns]=size(pros);

 if islogical(pros)
     flt = logical(zeros(9,1));%Mengecek image biner bukan dan memberinya 9 dimensi
 end
 
%Masukkan nilai-nilai dalam matriks 
for i=2:rows-1
    for j=2:columns-1
        flt=[pros(i-1,j-1),pros(i-1,j),pros(i-1,j+1),pros(i,j-1),pros(i,j),pros(i,j+1),pros(i+1,j-1),pros(i+1,j),pros(i+1,j+1)];
        pros(i,j)= median(flt);
    end
end

handles.gui.pros = pros;
handles.gui.I=handles.gui.pros;
imshow(handles.gui.I);
guidata(hObject,handles);

% --- Executes on button press in inverse.
function inverse_Callback(hObject, eventdata, handles)
% hObject    handle to inverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit_inverse = handles.gui.I;
[x, y] = size(edit_inverse);
if islogical (edit_inverse) %Mengecek image biner
    inverse=logical(zeros(x, y)); 
    for i=1:x 
        for j=1:y 
            if edit_inverse(i,j) ==0; %jika image bernilai 0 maka akan berubah menjadi 1
              inverse(i, j) = 1;
            else %sebaliknya
              inverse(i, j) = 0;
            end
        end
    end
elseif size(edit_inverse,3)==1 %Mengecek image geayscale
    inverse = uint8(zeros(x, y));
    for i=1:x 
        for j=1:y 
              inverse(i, j) = (255-edit_inverse(i, j)); %Melakukan aksi kebalikan 
        end
    end
elseif size(edit_inverse,3)==3 %Mengecek image RGB
    edit_inverse=uint8(edit_inverse);      
    inverse(:, :, 1)=uint8(bound(255-double(edit_inverse(:, :, 1)))); %R
    inverse(:, :, 2)=uint8(bound(255-double(edit_inverse(:, :, 2)))); %G
    inverse(:, :, 3)=uint8(bound(255-double(edit_inverse(:, :, 3)))); %B
    imshow(inverse);
end

handles.gui.inverse = inverse;
handles.gui.I = handles.gui.inverse;
imshow(handles.gui.I);
guidata(hObject,handles);


% --- Executes on button press in histogram.
function histogram_Callback(hObject, eventdata, handles)
% hObject    handle to histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.gui.I;

if size(I,3)==1; %jika image bertype grayscale
    [rows,cols] = size(I); 
    %Array histogram
    myhist = zeros(256,1);
    for k = 0:255
        myhist(k+1) = numel(find(I==k)); % jumlah elemen memiliki tingkat abu-abu sama dengan 'k'     
    end

    % menghitung cdf
    cdf = zeros(256,1);
    cdf(1) = myhist(1);
    for k = 2:256
        cdf(k) = cdf(k-1)+myhist(k);
    end

    % Mencari array histogram yang sama
    J = I;
    cumprob = cdf/(rows.*cols);
    equalizedhist = floor((cumprob).*255);

    for i = 1:rows
        for j = 1:cols
            for m = 0:255
                if (I(i,j)==m)
                    J(i,j) = equalizedhist(m+1);
                end
            end
        end
    end
handles.gui.J = J;
handles.gui.I = handles.gui.J;
imshow(handles.gui.I);
guidata(hObject, handles);
else    %jika image bertype RGB
    R = I(:,:,1);
    G = I(:,:,2);
    B = I(:,:,3);

    [rows,cols] = size(R);% Red

    % Array histogrm
    myhist = zeros(256,1);
    for k = 0:255
        myhist(k+1) = numel(find(R==k)); % jumlah elemen memiliki tingkat abu-abu sama dengan 'k'
    end

    % menghitung cdf
    cdf = zeros(256,1);
    cdf(1) = myhist(1);
    for k = 2:256
        cdf(k) = cdf(k-1)+myhist(k);
    end
    
    % Mencari array histogram yang sama
    R_histeq = R;
    cumprob = cdf/(rows.*cols);
    equalizedhist = floor((cumprob).*255);

    for i = 1:rows
        for j = 1:cols
            for m = 0:255
                if (R(i,j)==m)
                    R_histeq(i,j) = equalizedhist(m+1);
                end
            end
        end
    end


    % Mencari array histogram yang sama
    myeqhist = zeros(256,1);
    for k = 0:255
        myeqhist(k+1) = numel(find(R_histeq==k));
    end


    [rows,cols] = size(G); %Green

    % array Histogram 
    myhist = zeros(256,1);
    for k = 0:255
        myhist(k+1) = numel(find(G==k)); % jumlah elemen memiliki tingkat abu-abu sama dengan 'k'
    end

    % menghitung cdf
    cdf = zeros(256,1);
    cdf(1) = myhist(1);
    for k = 2:256
        cdf(k) = cdf(k-1)+myhist(k);
    end

    % Mencari array histogram yang sama
    G_histeq = G;
    cumprob = cdf/(rows.*cols);
    equalizedhist = floor((cumprob).*255);

    for i = 1:rows
        for j = 1:cols
            for m = 0:255
                if (G(i,j)==m)
                    G_histeq(i,j) = equalizedhist(m+1);
                end
            end
        end
    end


    % Mencari array histogram yang sama
    myeqhist = zeros(256,1);
    for k = 0:255
        myeqhist(k+1) = numel(find(G_histeq==k));
    end

    [rows,cols] = size(B);

    % array Histogram 
    myhist = zeros(256,1);
    for k = 0:255
        myhist(k+1) = numel(find(B==k));% jumlah elemen memiliki tingkat abu-abu sama dengan 'k'
    end

    % Menghitung cdf
    cdf = zeros(256,1);
    cdf(1) = myhist(1);
    for k = 2:256
        cdf(k) = cdf(k-1)+myhist(k);
    end

    % Mencari array histogram yang sama
    B_histeq = B;
    cumprob = cdf/(rows.*cols);
    equalizedhist = floor((cumprob).*255);

    for i = 1:rows
        for j = 1:cols
            for m = 0:255
                if (B(i,j)==m)
                    B_histeq(i,j) = equalizedhist(m+1);
                end
            end
        end
    end


    % Mencari array histogram yang sama
    myeqhist = zeros(256,1);
    for k = 0:255
        myeqhist(k+1) = numel(find(B_histeq==k));
    end
end
% Menyamkan gambar
RGB = cat(3,R_histeq,G_histeq,B_histeq);


handles.gui.RGB = RGB;
handles.gui.I = handles.gui.RGB;
imshow(handles.gui.I);
guidata(hObject, handles);

% --- Executes on button press in erosi.
function erosi_Callback(hObject, eventdata, handles)
% hObject    handle to erosi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit_erosi = handles.gui.I;
kernel =[1 1 1]; %karnel yang akan digunakan
C = padarray(edit_erosi,[0 1],1);%untuk mengganjal masukan citra secara explisit
erosi = false(size(edit_erosi));
for i=1:size(C,1)
    for j=1:size(C,2)-2
        R=C(i,j:j+2);
        K=find(kernel==1); %mencari karnel
       if(R(K)==1)
        erosi(i,j)=1;
        end
    end
end
handles.gui.erosi = erosi;
handles.gui.I = handles.gui.erosi;
imshow(handles.gui.I);
guidata(hObject, handles);

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I;
addpath(handles.gui.pathname);
I = imread(handles.gui.filename);
handles.gui.I = I;
imshow(I)
guidata(hObject,handles);

% --- Executes on button press in edge.
function edge_Callback(hObject, eventdata, handles)
% hObject    handle to edge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

edit_edge = handles.gui.I;
[rows,columns] =size(edit_edge);
if islogical(edit_edge) %jika gambar biner
    edit_edge1 = logical(zeros(rows,columns));
    edge = logical(zeros(8,1));%menggunakan 8 dimensi
    for i =1:rows
        for j =1:columns
           if j+1 > columns
               edge(1)=0;
           else
               edge(1) = edit_edge(i,j+1);
           end
            if i-1<=0||j+1>columns
                edge(2) =0;
            else 
                edge(2) = edit_edge(i-1,j+1);
            end
            if i-1<=0
                edge(3) = 0;
            else
                edge(3) = edit_edge(i-1,j);
            end
            if i-1<=0||j-1<=0
                edge(4) =0;
            else 
                edge(4) = edit_edge(i-1,j-1);
            end
            if j-1<=0
                edge(5) = 0;
            else
                edge(5) = edit_edge(i,j-1);
            end
            if i+1>rows||j-1<=0
                edge(6) =0;
            else 
                edge(6) = edit_edge(i+1,j-1);
            end
            if i+1>rows
                edge(7) = 0;
            else
                edge(7) = edit_edge(i+1,j);
            end
            if i+1>rows||j+1>columns
                edge(8) =0;
            else 
                edge(8) = edit_edge(i+1,j+1);
            end
            total = sum(edge);
            if total ==8
                edit_edge1(i,j)=0;
            else
                edit_edge1(i,j)=edit_edge(i,j);
            end
        end
    end
    handles.gui.edit_edge1 = edit_edge1;
    handles.gui.I=handles.gui.edit_edge1;
    imshow(handles.gui.I);
    guidata(hObject,handles);
else
    return;
end


% --- Executes on button press in histogram_btn.
function histogram_btn_Callback(hObject, eventdata, handles)
% hObject    handle to histogram_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = handles.gui.I;
if size(I,3)==1
    figure,
    subplot(1,1,1),imhist(I);
else
R =I(:,:,1); %Red
G =I(:,:,2); %Green
B =I(:,:,3); %Blue
    figure,
    subplot(1,3,1),histogram(R(:),256,'FaceColor','r','EdgeColor','r');
    subplot(1,3,2),histogram(G(:),256,'FaceColor','g','EdgeColor','g');
    subplot(1,3,3),histogram(B(:),256,'FaceColor','b','EdgeColor','b');
end
