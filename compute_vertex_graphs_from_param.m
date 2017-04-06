% param=load_tracking_param_from_pieces('X:\storage1\Antti\Exp186\0157_full_seq_boundary_correction\Exp186_full_F15_264_segmentation_parameters_');
% param.tracks=rmfield(param.tracks,'bounds');
% for c_ind=1:length(param.tracks)
%     param.tracks(c_ind).vertices=cell(size(param.tracks(c_ind).neighs));
%     for t_ind=1:length(param.tracks(c_ind).neighs)
%         param.tracks(c_ind).vertices{t_ind}=zeros(2,length(param.tracks(c_ind).neighs{t_ind}),'uint16');
%     end
%     param.tracks(c_ind).neighs=cell(size(param.tracks(c_ind).neighs));
%     for t_ind=1:length(param.tracks(c_ind).neighs)
%         param.tracks(c_ind).neighs{t_ind}=zeros(1,length(param.tracks(c_ind).neighs{t_ind}),'uint32');
%     end    
% end

for time=195:250
% for time=29:-1:3
% for time=11:200
% for time=100
%     disp(time)
    cents=zeros(length(param.tracks),4);
    c_ind=1;

    for i_ind=1:length(param.tracks)
        time_ind=find(param.tracks(i_ind).t==time);
        if isempty(time_ind)
            continue
        end

        cents(c_ind,:)=[double(param.tracks(i_ind).cent(time_ind,:)) i_ind time_ind];
        c_ind=c_ind+1;

    end

    cents(c_ind:end,:)=[];
    cents(end+1,:)=[0 0 0 0];



    %% input parameters
    tic
    % input image path
    s_img1=imread(['X:\storage1\Antti\Exp186\0157_full_seq_boundary_correction\segments\Exp186_full_F15_264_' num2str(time,'%.03d') '.tif']);

    % crop input image if needed

    if time<20
        inds_to_crop=[(500:1000:size(s_img1,2))-600;(500:1000:size(s_img1,2))+600];
    else
        inds_to_crop=[(250:500:size(s_img1,2))-300;(250:500:size(s_img1,2))+300];        
    end
    inds_to_crop(1)=1;
    inds_to_crop(end)=size(s_img1,2);
%     inds_to_crop=[1 ;size(s_img1,2)];

    for s_ind=1:size(inds_to_crop,2)
%     for s_ind=7
%         tic
        area_to_crop=[inds_to_crop(:,s_ind)'  1 size(s_img1,1)]; % [x x y y]
        s_img=s_img1(area_to_crop(3):area_to_crop(4),area_to_crop(1):area_to_crop(2));
        if s_ind>1
            left_bound=35;
        else
            left_bound=0;
        end
        if s_ind<size(inds_to_crop,2)
            right_bound=size(s_img,2)-35;
        else
            right_bound=size(s_img,2);
        end

        %% output parameters

        % centroids         : list of all centroids in segmentation [x_{c1}
        % y_{c1};x_{c2} y_{c2};...;x_{cN} y_{cN}], where c refers to centroid and N
        % is number of centroids. First centroid refers to region outside of
        % field of view.
        % vertex_coords     : list of all vertex points in segmentation [x_{v1}
        % y_{v1};x_{v2} y_{v2};...;x_{vM} y_{vM}], where v refers to vertex and M
        % is number of vertices
        % vertices          : list of vertices [v_{1,1} v_{1,2};...;v_{l,1}
        % v_{l,2}], where elements correspond to rows of the vertex_coords matrix,
        % each row represents an edge
        % neighs            : list of neighbouring centroids. Where elements
        % correspond to rows of the centroids matrix

        % 

        %%



        %% computation

        shrink=@(x)x(:);

        l_img1=zeros(size(s_img)+4);
        l_img1(3:end-2,3:end-2)=s_img;
        l_img1([2 end-1],2:end-1)=255;
        l_img1(2:end-1,[2 end-1])=255;
        
        l_img2=bwlabel(~l_img1,8);

        l_img1=bwmorph(l_img1,'skel',inf);

        % l_img1=l_img1-logical(conv2(double(l_img1),ones(2),'same')==4);

% six pixels together horizontal
        inds_to_go_through=find(conv2(double(l_img1),ones(2,3),'same')==6);        
        for i_ind=1:length(inds_to_go_through)

            if l_img1(inds_to_go_through(i_ind)-1)==0 & l_img1(inds_to_go_through(i_ind)+2)==1
                l_img1(inds_to_go_through(i_ind))=0;
            elseif l_img1(inds_to_go_through(i_ind)-1)==1 & l_img1(inds_to_go_through(i_ind)+2)==0
                l_img1(inds_to_go_through(i_ind)+1)=0;
            end
        end               

% six pixels together vertical
        inds_to_go_through=find(conv2(double(l_img1),ones(3,2),'same')==6);        
        for i_ind=1:length(inds_to_go_through)

            if l_img1(inds_to_go_through(i_ind)-size(l_img1,1))==0 & l_img1(inds_to_go_through(i_ind)+2*size(l_img1,1))==1
                l_img1(inds_to_go_through(i_ind))=0;
            elseif l_img1(inds_to_go_through(i_ind)-size(l_img1,1))==1 & l_img1(inds_to_go_through(i_ind)+2*size(l_img1,1))==0
                l_img1(inds_to_go_through(i_ind)+size(l_img1,1))=0;
            end
        end                       
        
% four pixels together vertical        
        inds_to_go_through=find(conv2(double(l_img1),ones(2),'same')==4);


       %%%%%%%


%         for i_ind=1:length(inds_to_go_through)
% 
%             if l_img1(inds_to_go_through(i_ind)-1-size(l_img1,1))==0
%                 l_img1(inds_to_go_through(i_ind))=0;
%             elseif l_img1(inds_to_go_through(i_ind)-1+2*size(l_img1,1))==0
%                 l_img1(inds_to_go_through(i_ind)+size(l_img1,1))=0;
%             elseif l_img1(inds_to_go_through(i_ind)+2+2*size(l_img1,1))==0
%                 l_img1(inds_to_go_through(i_ind)+1+size(l_img1,1))=0;
%             elseif l_img1(inds_to_go_through(i_ind)+2-size(l_img1,1))==0
%                 l_img1(inds_to_go_through(i_ind)+1)=0;
%             else
%                 l_img1(inds_to_go_through(i_ind))=0;
%                 if ~()
%                 l_img1(inds_to_go_through(i_ind)-1)=1;
%             end
%         end       
       
        for i_ind=1:length(inds_to_go_through)

            if l_img1(inds_to_go_through(i_ind)-1)~=0 | l_img1(inds_to_go_through(i_ind)-size(l_img1,1))~=0
                l_img1(inds_to_go_through(i_ind))=0;
            elseif l_img1(inds_to_go_through(i_ind)-1+size(l_img1,1))~=0 | l_img1(inds_to_go_through(i_ind)+2*size(l_img1,1))~=0
                l_img1(inds_to_go_through(i_ind)+size(l_img1,1))=0;
            elseif l_img1(inds_to_go_through(i_ind)+1+2*size(l_img1,1))~=0 | l_img1(inds_to_go_through(i_ind)+2+size(l_img1,1))~=0
                l_img1(inds_to_go_through(i_ind)+1+size(l_img1,1))=0;
            elseif l_img1(inds_to_go_through(i_ind)+2)~=0 | l_img1(inds_to_go_through(i_ind)+1-size(l_img1,1))~=0
                l_img1(inds_to_go_through(i_ind)+1)=0;
            else
                l_img1(inds_to_go_through(i_ind))=0;
                l_img1(inds_to_go_through(i_ind)-1)=1;
            end
        end       

       %%%%%%%
        
        
        
                                    l_img1=bwlabel(~l_img1,4);

                                    match_pairs=unique([l_img2(:) l_img1(:)],'rows');
                                    subs_a=match_pairs(match_pairs(:,1)>0 & match_pairs(:,2)>0 & match_pairs(:,1)~=match_pairs(:,2),:);
                                    subs_a=sortrows(subs_a,[2 1]);

                                    l_img3=l_img1;
                                    for i_ind=1:size(subs_a,1)
                                        l_img3(l_img1(:)==subs_a(i_ind,2))=subs_a(i_ind,1);
                                    end
                                    mapping_subs=[unique(l_img3(:)) (0:length(unique(l_img3(:)))-1)'];
                                    mapping_subs(mapping_subs(:,1)==mapping_subs(:,2),:)=[];
                                    for i_ind=1:size(mapping_subs,1)
                                        l_img3(l_img3(:)==mapping_subs(i_ind,1))=mapping_subs(i_ind,2);
                                    end        

                                    l_img1=l_img3;
        

        num_segments=length(unique(l_img1(:)))-1; % first segment is outside of the region

        centroids=regionprops(l_img1,'centroid');
        centroids=reshape(struct2array(centroids)',2,length(centroids))';
        centroids=centroids-2;
        centroids(1,:)=nan;

        [vals, cent_inds]=min((repmat(centroids(:,1),[1 size(cents,1)])-repmat(cents(:,1)'+1-area_to_crop(1),[size(centroids,1) 1])).^2+(repmat(centroids(:,2),[1 size(cents,1)])-repmat(cents(:,2)'+1-area_to_crop(3),[size(centroids,1) 1])).^2);
        [vals1, cent_inds1]=min((repmat(centroids(:,1),[1 size(cents,1)])-repmat(cents(:,1)'+1-area_to_crop(1),[size(centroids,1) 1])).^2+(repmat(centroids(:,2),[1 size(cents,1)])-repmat(cents(:,2)'+1-area_to_crop(3),[size(centroids,1) 1])).^2,[],2);
        
        cent_inds1(1)=size(cents,1); % boundary indicator in tracklet structure

        % close all
        % figure(1)
        % clf
        % imagesc(l_img1)

        b_label=l_img1(1,1);

        b_inds=find(l_img1(:)==0);

        b_coords=[floor((b_inds-1)/size(l_img1,1))+1 mod(b_inds-1,size(l_img1,1))+1];

        eigth_pixels=[l_img1(b_inds-size(l_img1,1)-1) l_img1(b_inds-size(l_img1,1)) l_img1(b_inds-size(l_img1,1)+1) ...
            l_img1(b_inds-1) l_img1(b_inds+1) ...
            l_img1(b_inds+size(l_img1,1)-1) l_img1(b_inds+size(l_img1,1)) l_img1(b_inds+size(l_img1,1)+1)];

        eigth_pixels=sort(eigth_pixels,2);

        eigth_pixels_diffs=diff(eigth_pixels,1,2)>0;

        one_cell=sum(eigth_pixels_diffs,2)==1;
        one_cell_diff=eigth_pixels_diffs(one_cell,:);
        one_cell_diff=[one_cell_diff ones(size(one_cell_diff,1),1)]';
        one_cell_actual=eigth_pixels(one_cell,:)';
        one_cell_inds=reshape(one_cell_actual(logical(one_cell_diff(:))),[2 size(one_cell_diff,2)])';
        one_cell_coords=b_coords(one_cell,:);

        two_cells=sum(eigth_pixels_diffs,2)==2;
        two_cells_diff=eigth_pixels_diffs(two_cells,:);
        two_cells_diff=[two_cells_diff ones(size(two_cells_diff,1),1)]';
        two_cells_actual=eigth_pixels(two_cells,:)';
        two_cells_inds=reshape(two_cells_actual(logical(two_cells_diff(:))),[3 size(two_cells_diff,2)])';
        two_cells_coords=b_coords(two_cells,:);

        three_cells=sum(eigth_pixels_diffs,2)==3;
        three_cells_diff=eigth_pixels_diffs(three_cells,:);
        three_cells_diff=[three_cells_diff ones(size(three_cells_diff,1),1)]';
        three_cells_actual=eigth_pixels(three_cells,:)';
        three_cells_inds=reshape(three_cells_actual(logical(three_cells_diff(:))),[4 size(three_cells_diff,2)])';
        three_cells_coords=b_coords(three_cells,:);

        four_cells=sum(eigth_pixels_diffs,2)==4;
        four_cells_diff=eigth_pixels_diffs(four_cells,:);
        four_cells_diff=[four_cells_diff ones(size(four_cells_diff,1),1)]';
        four_cells_actual=eigth_pixels(four_cells,:)';
        four_cells_inds=reshape(four_cells_actual(logical(four_cells_diff(:))),[5 size(four_cells_diff,2)])';
        four_cells_coords=b_coords(four_cells,:);

        vertex_coords=[three_cells_coords ;four_cells_coords]; %[x y ] [three neighbours, four neighbours]

        % if size(one_cell_inds,1)+size(two_cells_inds,1)+size(three_cells_inds,1)+size(four_cells_inds,1)~=size(eigth_pixels,1)
        %     disp('number of neighs not matching')
        % end
        % if size(one_cell_inds,1)>0
        %    disp('one neighbours found') 
        % end
        % if size(four_cells_inds,1)>0
        %    disp('four neighbours found') 
        % end

        vertices=zeros(size(vertex_coords,1)*3,2);
        vertices_ind=1;
        neighs=zeros(size(vertex_coords,1)*3,2);
        neighs_ind=1;
        vertex_coords=vertex_coords-2;

        vertex_coords(:,1)=vertex_coords(:,1)+area_to_crop(1)-1;
        vertex_coords(:,2)=vertex_coords(:,2)+area_to_crop(3)-1;
        vertex_coords(end+1,:)=[inf inf];

        for c_ind=2:num_segments
        % for c_ind=1:23

            cell_two_inds=two_cells_inds(two_cells_inds(:,2)==c_ind |two_cells_inds(:,3)==c_ind,:);
            cell_two_inds(cell_two_inds==c_ind)=0;
            cell_two_inds=sort(cell_two_inds,2);
            cell_two_inds=sortrows([two_cells_coords(two_cells_inds(:,2)==c_ind |two_cells_inds(:,3)==c_ind,:) cell_two_inds(:,3)],3);

            cell_three_inds=three_cells_inds(sum(three_cells_inds==c_ind,2)>0,:);
            cell_three_inds(cell_three_inds==c_ind)=0;
            cell_three_inds=sort(cell_three_inds,2);
            cell_three_inds=[cell_three_inds find(sum(three_cells_inds==c_ind,2)>0)];
            cell_three_inds=sortrows([three_cells_coords(sum(three_cells_inds==c_ind,2)>0,:) cell_three_inds(:,3:5)],[3 4 5]);

            cell_four_inds=four_cells_inds(sum(four_cells_inds==c_ind,2)>0,:);
            cell_four_inds(cell_four_inds==c_ind)=0;
            cell_four_inds=sort(cell_four_inds,2);
            cell_four_inds=[cell_four_inds find(sum(four_cells_inds==c_ind,2)>0)];
            if ~isempty(cell_four_inds)
                cell_four_inds=sortrows([four_cells_coords(sum(four_cells_inds==c_ind,2)>0,:) cell_four_inds(:,3:6)],[3 4 5 6]);
            else
                cell_four_inds=zeros(0,6);
            end

            cell_three_four_inds=[[cell_three_inds(:,1:4) zeros(size(cell_three_inds,1),1) cell_three_inds(:,5)];cell_four_inds];
            cell_three_four_inds((size(cell_three_inds,1)+1):end,6)=cell_three_four_inds((size(cell_three_inds,1)+1):end,6)+size(three_cells_inds,1);

            all_coords=[[cell_two_inds(:,1:2) 2*ones(size(cell_two_inds,1),1) (1:size(cell_two_inds,1))'];...
                [cell_three_four_inds(:,1:2) 34*ones(size(cell_three_four_inds,1),1) (1:size(cell_three_four_inds,1))']]; %[x y num_points ind_to_struct]


            cur_point=1;
            cur_coord=all_coords(cur_point,:);

            if isempty(cell_three_four_inds)
                lines_to_add=[ 0 0 cell_two_inds(cur_coord(4),3) 0 0 size(vertex_coords,1) 0 0 cell_two_inds(cur_coord(4),3) 0 0 size(vertex_coords,1)];
                lines_to_add_check=1;
                cur_end0=cur_coord;
            else
                lines_to_add=[];  
                lines_to_add_check=[];
                cur_end0=all_coords(size(cell_two_inds,1)+1,:);
            end

            n=zeros(10,12);
            n_ind=1;

        %     centroids(1,:)=nan;

            three_four_ind_found=0;
            while 1
%                 plot(all_coords(cur_point,1),all_coords(cur_point,2),'.g');
                all_coords(cur_point,:)=[];

                if size(all_coords,1)==0
                    break;
                end

                diffs=((all_coords(:,1)-cur_coord(1)).^2+(all_coords(:,2)-cur_coord(2)).^2);
                
% %                 inds_1 = find(diffs==1);
% %                 if size(inds_1,1)>0
% %                     d=1;
% %                     cur_point=inds_1(1);
% %                 else
% %                     inds_2 = find(diffs==2);
% %                     if size(inds_2,1)==1
% %                         d=2;
% %                         cur_point=inds_2(1);
% %                     elseif size(inds_2,1)>1
% %                         for ind_ind=1:size(inds_2,1)
% %                             if l_img2(all_coords(inds_2(ind_ind),2),cur_coord(1))==0 | l_img2(cur_coord(2),all_coords(inds_2(ind_ind),1))==0
% %                                 d=2;
% %                                 cur_point=inds_2(ind_ind);
% %                                 break;
% %                             end
% %                         end
% %                     else
% %                         [d,cur_point]=min((all_coords(:,1)-cur_coord(1)).^2+(all_coords(:,2)-cur_coord(2)).^2);
% %                     end
% %                 end
                

                inds_1 = find(diffs==1);
                inds_2 = find(diffs==2);

% choose shortest path if only one neighbour was found                
                if size(inds_1,1)==1 & size(inds_2,1)==0
                    d=1;
                    cur_point=inds_1(1);
% choose shortest path if only one neighbour was found                    
                elseif size(inds_1,1)==0 & size(inds_2,1)==1
                    d=2;
                    cur_point=inds_2(1);
% choose longer path if no neighbours was found
                elseif size(inds_1,1)==0 & size(inds_2,1)==0
                    [d,cur_point]=min((all_coords(:,1)-cur_coord(1)).^2+(all_coords(:,2)-cur_coord(2)).^2);
% otherwise choose neighbouring pixel with least neighbours
                else
                    inds_12=[inds_1 ones(size(inds_1));inds_2 ones(size(inds_2))*2];
                    inds_12(inds_12(:,2)==1,4)=1;
                    for ind_ind=1:size(inds_12,1)
                        inds_12(ind_ind,3)=sum(((all_coords(:,1)-all_coords(inds_12(ind_ind,1),1)).^2+(all_coords(:,2)-all_coords(inds_12(ind_ind,1),2)).^2)<=2);
                        if inds_12(ind_ind,2)==2 
                            inds_12(ind_ind,4)=l_img2(all_coords(inds_12(ind_ind,1),2),cur_coord(1))==0 | l_img2(cur_coord(2),all_coords(inds_12(ind_ind,1),1))==0;
                        end
                    end
                    
                    sorted_inds_12=sortrows(inds_12,[2 3]);
                    ind_to_use=find(sorted_inds_12(:,3)==2 & sorted_inds_12(:,4)==1,1,'first');
                    if isempty(ind_to_use)
                        ind_to_use=find(sorted_inds_12(:,4),1,'first');
                    end
                    d=sorted_inds_12(ind_to_use,2);
                    cur_point=sorted_inds_12(ind_to_use,1);
                    
                end
                  
                %%%

                
                cur_coord=all_coords(cur_point,:);
                
%                 plot(cur_coord(1),cur_coord(2),'.g')
%                 pause
                
                if d>2 % a hole inside of a segment
    %                 disp('a hole inside of a segment')                
%                     break

                    lines_to_add(end+1,:)=[ 0 0 cell_two_inds(cur_coord(4),3) 0 0 size(vertex_coords,1) 0 0 cell_two_inds(cur_coord(4),3) 0 0 size(vertex_coords,1)];
                    
                    lines_to_add_check(end+1)=three_four_ind_found;
                    three_four_ind_found=0;
                    
                    % % %                     cur_end0=cur_coord;
                elseif cur_coord(3)==34   
                    n(n_ind,:)=[ cell_three_four_inds(cur_end0(4),:) cell_three_four_inds(cur_coord(4),:)];
                    n_ind=n_ind+1;
                    cur_end0=cur_coord;
                    three_four_ind_found=1;
%                     pause
                end
            end

            if n_ind>1
                n(1,1:6)= n(n_ind-1,7:12);
            end

            n((n_ind):end,:)=[];

            if ~isempty(lines_to_add)
                lines_to_add(logical([lines_to_add_check(2:end) three_four_ind_found]),:)=[];
            end
            
            for l_ind=size(lines_to_add,1):-1:1
                if sum(sum(n(:,[3 4 5])==lines_to_add(l_ind,3)))>0
                    lines_to_add(l_ind,:)=[];
                end
            end
                
            n=[n; lines_to_add];
            

            n(n(:,3)==n(:,9)& n(:,3)>0,13)=n(n(:,3)==n(:,9) & n(:,3)>0,3);
            n(n(:,3)==n(:,10)& n(:,3)>0,13)=n(n(:,3)==n(:,10) & n(:,3)>0,3);
            n(n(:,3)==n(:,11)& n(:,3)>0,13)=n(n(:,3)==n(:,11) & n(:,3)>0,3);
            n(n(:,4)==n(:,9)& n(:,4)>0,13)=n(n(:,4)==n(:,9) & n(:,4)>0,4);
            n(n(:,4)==n(:,10)& n(:,4)>0,13)=n(n(:,4)==n(:,10) & n(:,4)>0,4);
            n(n(:,4)==n(:,11)& n(:,4)>0,13)=n(n(:,4)==n(:,11) & n(:,4)>0,4);
            n(n(:,5)==n(:,9)& n(:,5)>0,13)=n(n(:,5)==n(:,9) & n(:,5)>0,5);
            n(n(:,5)==n(:,10)& n(:,5)>0,13)=n(n(:,5)==n(:,10) & n(:,5)>0,5);
            n(n(:,5)==n(:,11)& n(:,5)>0,13)=n(n(:,5)==n(:,11) & n(:,5)>0,5);

            % plot(n(:,[1 6])',n(:,[2 7])','-r')

            vertices(vertices_ind:(vertices_ind+size(n,1)-1),:) = n(:,[6 12]); % indices to vertex_coords
            neighs(neighs_ind:(neighs_ind+size(n,1)-1),:) = [ones(size(n,1),1)*c_ind n(:,13)]; % indices to centroids

            vertices_ind=vertices_ind+size(n,1);
            neighs_ind=neighs_ind+size(n,1);

%             ind=find(cent_inds==c_ind & vals<=1);
            ind=find(cent_inds==c_ind & centroids(c_ind,1)>=left_bound & centroids(c_ind,1)<=right_bound);

            if size(ind,2)==0
                continue
            end
            
            if size(ind,2)>1                
                [~,m_ind]=min(vals(ind));   
                ind=ind(m_ind);             
                [time s_ind c_ind]          
%                 plot(centroids(c_ind,1),centroids(c_ind,2),'ow')
%                 continue
            end                             
            
%             left_bound
%             right_bound

            cell_i=cents(ind,3);
            cell_time=cents(ind,4);

            param.tracks(cell_i).neighs{cell_time}=uint32(cents(cent_inds1(n(:,13)),3)');    
            param.tracks(cell_i).vertices{cell_time}=uint16(vertex_coords(n(:,6),:)');


        %     close all
        %     figure(1)
        %     hold on
        %     clf
        %     imagesc(l_img1)

        %     hold on
        % 
        %     plot([centroids(c_ind*ones(size(neighs,1),1),1) centroids(neighs,1)]',[centroids(c_ind*ones(size(neighs,1),1),2) centroids(neighs,2)]','+-k','linewidth',2)
        % 
        %     plot([vertex_coords(vertices(:,1),1) vertex_coords(vertices(:,2),1)]',[vertex_coords(vertices(:,1),2) vertex_coords(vertices(:,2),2)]','o-w','linewidth',2)

        end

        vertices(vertices_ind:end,:)=[];
        neighs(neighs_ind:end,:)=[];

        temp=unique([sort(vertices,2) sort(neighs,2)],'rows');
        vertices=temp(:,1:2);
        neighs=temp(:,3:4);


%         disp([s_ind toc])

    end

    for c_ind=1:length(param.tracks)
        t_ind=param.tracks(c_ind).t==time;
        if sum(t_ind)==0
            continue
        end
        if sum(sum(param.tracks(c_ind).vertices{t_ind}==0)>1)>0
            disp([c_ind find(t_ind)])
        end
    end
    
    disp([time toc])
    pause(2)
    
    
    
end


%% save data

%% saving            
file_name=[param.seq_name '_segmentation_parameters_vertices_1_250_'];
% file_name=[param.seq_name '_segmentation_parameters_vertices_1_194_'];

param1=[];   
file_num=1;

track_interval=[1 length(param.tracks); length(param.tracks) length(param.tracks)];
while 1                    
    param1.tracks=param.tracks(track_interval(file_num,1):track_interval(file_num,2));
    size_param=whos('param1');
    size_param.bytes
    while size_param.bytes>(2e9)
        track_interval(file_num,2)=round(track_interval(file_num,2)*0.9);
        track_interval(file_num+1,1)=track_interval(file_num,2)+1;

        param1.tracks=param.tracks(track_interval(file_num,1):track_interval(file_num,2));
        size_param=whos('param1');                            
    end

    param1.track_interval=track_interval(file_num,:);
    save([file_name num2str(file_num)],'param1');                       

    if track_interval(file_num,2)==length(param.tracks)
        break;
    end

    file_num=file_num+1;
    track_interval(end+1,:)=[length(param.tracks) length(param.tracks)];
end

clear('param1');

f_names=fieldnames(param);
f_names(strcmp(f_names,'tracks'))=[];        
param0=[];
for i_ind=1:length(f_names)
    f_name=f_names(i_ind);
    f_name=f_name{1};        
    param0=setfield(param0,f_name,getfield(param,f_name));
end    
save([file_name num2str(0)],'param0');   

clear('param0'); 

            
%% end save data

% plot(b_coords(two_cells,1),b_coords(two_cells,2),'.r');
% 
% plot(b_coords(three_cells,1),b_coords(three_cells,2),'or');

%% visualisation (this can be commented out)

% % [toc numel(s_img)/1e6]
% % 
% % 
% % figure
% % imagesc(s_img)
% % hold on
% % 
% % plot([centroids(neighs(:,1),1) centroids(neighs(:,2),1)]',[centroids(neighs(:,1),2) centroids(neighs(:,2),2)]','+-g','linewidth',2)
% % 
% % plot([vertex_coords(vertices(:,1),1) vertex_coords(vertices(:,2),1)]',[vertex_coords(vertices(:,1),2) vertex_coords(vertices(:,2),2)]','o-w','linewidth',2)
% % 
% % axis equal
% % 
% % 
% % %%
% % 
% % 
% % %% 
% % cell_time=1;
% % 
% % cents=zeros(length(param.tracks),3);
% % c_ind=1;
% % 
% % for i_ind=1:length(param.tracks)
% %     if param.tracks(i_ind).t(1)==1 & param.tracks(i_ind).cent(1,1)>area_to_crop(1) & param.tracks(i_ind).cent(1,2)>area_to_crop(3) & param.tracks(i_ind).cent(1,1)<area_to_crop(2) & param.tracks(i_ind).cent(1,2)<area_to_crop(4)
% %         cents(c_ind,:)=[param.tracks(i_ind).cent(1,:) i_ind];
% %         c_ind=c_ind+1;
% %     end
% % end
% % 
% % cents(c_ind:end,:)=[];
% % 
% % [vals, cent_inds]=min((repmat(centroids(:,1),[1 size(cents,1)])-repmat(cents(:,1)'+1-area_to_crop(1),[size(centroids,1) 1])).^2+(repmat(centroids(:,2),[1 size(cents,1)])-repmat(cents(:,2)'+1-area_to_crop(3),[size(centroids,1) 1])).^2);
% % 
% % for c_ind=1:size(cents,1)
% %     n_inds_to_consider=find(logical(sum(neighs==cent_inds(c_ind),2)));
% %     
% %     for n_ind=n_inds_to_consider';
% %         cents_temp=centroids(neighs(n_ind,:),:);
% %         if sum(isnan(cents_temp(:)))>0
% %             continue
% %         end        
% %         
% %         cell_id1=cents(find(cent_inds==(neighs(n_ind,logical(1-(neighs(n_ind,:)==cent_inds(c_ind)))))),3);
% %         if isempty(cell_id1)
% %             continue
% %         end
% %         cell_id2=cents(c_ind,3);
% %         
% %         if sum(param.tracks(cell_id1).neighs{cell_time}==cell_id2)==0
% %             disp([cell_id1 cell_id2])
% %             continue
% %         end
% %         
% %         param.tracks(cell_id1).vertices{cell_time}(1:4,param.tracks(cell_id1).neighs{cell_time}==cell_id2)=uint16([vertex_coords(vertices(n_ind,1),:)'-1+area_to_crop([1 3])' ; vertex_coords(vertices(n_ind,2),:)'-1+area_to_crop([1 3])']); % [x y x y]'
% %         param.tracks(cell_id2).vertices{cell_time}(1:4,param.tracks(cell_id2).neighs{cell_time}==cell_id1)=uint16([vertex_coords(vertices(n_ind,1),:)'-1+area_to_crop([1 3])' ; vertex_coords(vertices(n_ind,2),:)'-1+area_to_crop([1 3])']); % [x y x y]'
% %         
% %     end    
% % end













