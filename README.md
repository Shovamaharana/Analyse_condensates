# nuclear_condensates
This is a Matlab Code to analyse nuclear condensate images for enrichment of respective markers. We start here with the best focussed images of the condesates and the respective images in the dapi channel. We keep all the channel files in one folder. The channel names are usually separated by '_c001.tif', '_c002.tif' or 'c:1-3.tif', 'c:2-3.tif' etc.
Each filed of view opens one after another. We use the dapi channel to automatically detect the nucleus. We separate the nucleus and look for intranuclear condensates.
We use a 12 pixel rim around the nucleus to calculate cytoplasmic intensity.
The results are compiled in variable-"all_result1", which can be copied into an excel.


clc;
clear all;

%% make directory of the files needs to be analysed
dv_folder = '/Volumes/Backup Plus/shrna_fibrillarin';
cd(dv_folder);

file = dir('c:2-3*');
N=size(file);

%% loop for all files
all_result=[];all_result1=[];result_indipml=[];

for a=1:N
    
    file(a).name;
    nm_C2 = uint16(imread(file(a).name));a
    
    
    str1 = file(a).name;
    str2 = regexp (str1,'c:2-3','split'); % nucleoulus/pml/cajal channel
    str = str2{1,2};
    
    str02 = strcat('c:1-3',str);%%% dapi channel
    nm_C1 = uint16(imread(str02));
  

    %% loop for each cell in the image
     all_cyto_mask=[];all_nuc_mask=[];sum_cytomask = uint16(ones(size(nm_C1,1),size(nm_C1,2)));
  
   blue = uint16(zeros(size(nm_C1,1),size(nm_C1,2)));
   RGB= cat(3,(nm_C2),blue,blue);
   
    img_dapi = nm_C1;
    se = strel('disk',15)
    background = imopen(img_dapi,se);
    img_dapi2 = img_dapi - background;
    imshow(img_dapi);
    
    img_dapi3 = imadjust(img_dapi2);
    
 
       
        
        % create a freehand mask to select a few nuclei to set threshold for auto sepration of nuclei, preferably select a nuclei with less Dapi intensity
        RGB =RGB.*sum_cytomask;
        img=imagesc(nm_C1);title('Select cyto Boundary');
        %img=imshow(RGB);title('Select cyto Boundary');
        h = imfreehand;

        cyto_mask0 = createMask(h, img);
        invcyto_mask = uint16(imcomplement(cyto_mask0));
        cyto_mask = uint16(cyto_mask0);
        close all; figure (1);imagesc(cyto_mask);pause(0.1)
        all_cyto_mask=cat(3,all_cyto_mask,cyto_mask0);
         % compute for flag
        lin_cell=double(cyto_mask(cyto_mask>0));
        sum_area=sum(lin_cell);
       
            
            % nucleus mask by thresholding DAPI
            img_cell=uint16(cyto_mask.*nm_C1);
            lin_I1cell=double(img_cell(img_cell>0));
            mn1 = mean (lin_I1cell);
            sd1 = std(lin_I1cell);
            
            Iblur1 = imgaussfilt(nm_C1,4);
            
            nuc_thr = Iblur1 > ((0.65 *mn1)+(0.9*sd1));%%%% select nucleus
            nuc_thr = bwareaopen(nuc_thr,5000);
            SE = strel("disk",4)
            nuc_thr = imdilate(nuc_thr, SE);
            nuc_thr = imerode(nuc_thr, SE);
            nuc_thr = imfill(nuc_thr,'holes');
            imagesc(nuc_thr);pause(0.3)
            % detect individual nuclei
            cc = bwconncomp(nuc_thr,4);
            
            cc.NumObjects
            %%
            %labeled = labelmatrix(cc);
            %RGB_label = label2rgb(labeled,'spring','c','shuffle');
            %imshow(RGB_label)
            
            for b =1:cc.NumObjects
                grain = false(size(nuc_thr));
                grain(cc.PixelIdxList{b}) = true;
                nuc_mask =uint16(grain);
            %imagesc(grain);pause(0.1);
            
            % calculate cytoplasm mask by dilating nuclear mask
           SE = strel("disk",12);
           nuc_expand = imdilate(grain,SE);
           nuc_rim = nuc_expand - grain;
           imagesc(nuc_rim);
           cytoplasm=uint16(nuc_rim);
            
             % speckle or nucleouls mask by thresholding nuclear marker %%%% 
             img_cell=uint16(nuc_mask.*nm_C2);
             lin_I1cell=double(img_cell(img_cell>0));
             mn2 = mean (lin_I1cell);
             sd2 = std(lin_I1cell);
             Iblur1 = imgaussfilt(img_cell,2);
            sp_thr = Iblur1 > ((1.12*mn2)+(0.15*sd2));
    
            sp_thr1 = bwareaopen(sp_thr,30);
            all_sp_mask = uint16(sp_thr1); %%%% select speckles
            %figure (2);imagesc( all_sp_mask);title('speckle');pause(0.3);
            % cytoplasm with no foci and no sg
            sp_thri =uint16(imcomplement(sp_thr1));
            nuc_nosp_mask = nuc_mask .*  sp_thri;
            %figure (1);imagesc(nuc_nosp_mask);title('no sp nuc');pause(0.3);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
            nucleoplasm = nuc_nosp_mask ; %%% nucleoplasm
             nucleopla_int1 = nucleoplasm.*nm_C2; %%%% real image quanti
                 lin_I1nucleoplasm1=double(nucleopla_int1(nucleopla_int1>0));
                   mn2_nucleoplasm = mean (lin_I1nucleoplasm1);
                   sd2_nucleoplasm = std(lin_I1nucleoplasm1);
                   nuc_area = sum(sum(nuc_mask));
                    
           cytopla2_int =cytoplasm .*nm_C2; %%%% real image quati J2
                   lin_I2cytoplasm=double(cytopla2_int(cytopla2_int>0));
                   mn2_cytoplasm = mean (lin_I2cytoplasm);
                   sd2_cytoplasm = std(lin_I2cytoplasm);
 
                   sp_int = all_sp_mask .*nm_C2;
                   lin_I2sp =double(sp_int(sp_int > 0));
                   mn2_sp = mean(lin_I2sp); %%%%nucleolus
                   sd2_sp = std(lin_I2sp);
                   sp_area = sum(sum(all_sp_mask));
   
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
           result1=cat(2,mn2_nucleoplasm,sd2_nucleoplasm, nuc_area,mn2_cytoplasm,sd2_cytoplasm,mn2_sp,sd2_sp, sp_area);
          all_result1=cat(1,all_result1,result1);
       
           close all;
                
            end
                
            
            
       
   
   
   
end
    
