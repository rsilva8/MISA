% load('MISA_MISA_ERPfMRI.mat')

%% ERP Components correspond 1-to-1 (GT vs Est.)
% Ranges of -.1 to .21 in y-axis suffice
figure('position',[131   750   560   420/2])
plot(A{1}(:,1),'b','LineWidth',2)
hold on
plot(Wout{1}(1,:),':r','LineWidth',2)
axis tight
ylim([-.1 .21])
plot_crosshair(gca)
pos_ = get(gca,'position');
pos_(4) = 2*0.34116; set(gca,'position', pos_);

figure('position',[131   450   560   420/2])
plot(A{1}(:,2),'b','LineWidth',2)
hold on
plot(Wout{1}(2,:),':r','LineWidth',2)
axis tight
ylim([-.1 .21])
plot_crosshair(gca)
pos_ = get(gca,'position');
pos_(4) = 2*0.34116; set(gca,'position', pos_);

[u, s, v] = svd(Wout{1}(3:4,:),'econ');
figure('position',[131   50   560   420])
subplot(211)
plot(A{1}(:,3),'b','LineWidth',2)
hold on
plot(Wout{1}(3,:),':r','LineWidth',2)
plot(-v(:,1)',':g','LineWidth',2)
axis tight
ylim([-.1 .21])
plot_crosshair(gca)
subplot(212)
plot(A{1}(:,4),'b','LineWidth',2)
hold on
plot(Wout{1}(4,:),':r','LineWidth',2)
plot(v(:,2)',':g','LineWidth',2)
axis tight
ylim([-.1 .21])
plot_crosshair(gca)

%% Need to load mask from here:
% cd('\\loki\export\mialab\users\rsilva\projects\MultivariateICA\MISA\fresh\data\Wgt\fMRI')
% msk=spm_read_vols(spm_vol('rest_hcpMask.img'));
% load('\\loki\export\apps\linux-x86\matlab\toolboxes\MCIv4\CM_coldhot_256.mat', 'CM')

for cc = 1:6
    b = zeros(size(msk));
    b(logical(msk)) = A{2}(:,cc);
    [ii(cc), jj(cc), kk(cc)] = ...
        ind2sub(size(b),find( abs(b) == max(abs(b(:))) ));
end
ii(1) = 21; jj(1) = 50; kk(1) = 33;
% ii(2) = ii(3); jj(2) = jj(3); kk(2) = kk(3);
% ii(6) = ii(5); jj(6) = jj(5); kk(6) = kk(5);
ii(6) = 15; jj(6) = 13; kk(6) = 5;

bix = [1 2 3 4 5 6];
% Wout{2}(1,:) = -Wout{2}(1,:);
Wout{2}(2,:) = -Wout{2}(2,:);
Wout{2}(5,:) = Wout{2}(6,:);
Wout{2}(6,:) = Wout{2}(5,:);
[u1, s1, v1_] = svd(Wout{2}(2:3,:),'econ');
[W1,wht]=icatb_runica(v1_', 'sphering', 'none', 'verbose', 'off');
v1 = W1*v1_';
v1 = max(abs(b(:)))*diag(1./max(abs(v1),[],2))*v1;
[~,ix] = max(reshape(corr(A{2}(:,2:3),v1'),[],1));
[I,J] = ind2sub([2 2],ix);
if I ~= J
    v1(1:end,:) = v1(end:-1:1,:);
end
% [u2, s2, v2_] = svd(Wout{2}(5:6,:),'econ');
% [W2,wht]=icatb_runica(v2_', 'sphering', 'none', 'verbose', 'off');
% v2 = W2*v2_';
v2 = Wout{2}(5:6,:);

no3rdimg = [1 4 5 6];

the_vol = spm_vol('rest_hcpMask.img');

for cc = 1:6
%     for dd = 1:3
        a = zeros(size(msk));
        b = zeros(size(msk));
        c = zeros(size(msk));
        a(logical(msk)) = Wout{2}(bix(cc),:);
        b(logical(msk)) = A{2}(:,cc);
        if cc == 2 || cc == 3
            c(logical(msk)) = v1(cc-1,:);
        elseif cc == 5 || cc == 6
            c(logical(msk)) = v2(cc-4,:);
        end
        the_vol.fname = ['MISAfusion_fMRI' num2str(cc) '_GT.nii'];
        the_vol.private.dat.fname = the_vol.fname;
        spm_write_vol(the_vol,b);
        the_vol.fname = ['MISAfusion_fMRI' num2str(cc) '_est.nii'];
        the_vol.private.dat.fname = the_vol.fname;
        if sum((no3rdimg - cc) == 0) > 0
            spm_write_vol(the_vol,a);
        else
            spm_write_vol(the_vol,c);
        end
%         switch dd
%             case 1
%                 if sum((no3rdimg - cc) == 0) > 0
%                     A_ = [rot90(squeeze(b(ii(cc),:,:))',2)...
%                         rot90(squeeze(a(ii(cc),:,:))',2)];
%                 else
%                     A_ = [rot90(squeeze(b(ii(cc),:,:))',2)...
%                         rot90(squeeze(a(ii(cc),:,:))',2) ...
%                         rot90(squeeze(c(ii(cc),:,:))',2)];
%                 end
%             case 2
%                 if sum((no3rdimg - cc) == 0) > 0
%                     A_ = [rot90(squeeze(b(:,jj(cc),:))',2)...
%                         rot90(squeeze(a(:,jj(cc),:))',2)];
%                 else
%                     A_ = [rot90(squeeze(b(:,jj(cc),:))',2)...
%                         rot90(squeeze(a(:,jj(cc),:))',2) ...
%                         rot90(squeeze(c(:,jj(cc),:))',2)];
%                 end
%             case 3
%                 if sum((no3rdimg - cc) == 0) > 0
%                     A_ = [rot90(squeeze(b(:,:,kk(cc)))',2)...
%                         rot90(squeeze(a(:,:,kk(cc)))',2)];
%                 else
%                     A_ = [rot90(squeeze(b(:,:,kk(cc)))',2)...
%                         rot90(squeeze(a(:,:,kk(cc)))',2) ...
%                         rot90(squeeze(c(:,:,kk(cc)))',2)];
%                 end
%         end
%         figure
%         imagesc(A_,max(abs(b(:)))*[-1 1])
%         colormap(CM)
%         axis equal tight
%         drawnow()
%     end
end







