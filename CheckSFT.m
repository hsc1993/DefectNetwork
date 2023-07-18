function check=CheckSFT(rn,rn0,doSFT,SFT_plane)

SFTol=-1e-4;
for i=0:(size(SFT_plane,1)/4)-1
    SFT_index=4*i+1;
    if((dot((rn0-(SFT_plane(SFT_index,4:6))),SFT_plane(SFT_index,1:3))<SFTol)&&...
      (dot((rn0-(SFT_plane(SFT_index+1,4:6))),SFT_plane(SFT_index+1,1:3))<SFTol)&&...        
      (dot((rn0-(SFT_plane(SFT_index+2,4:6))),SFT_plane(SFT_index+2,1:3))<SFTol)&&...
      (dot((rn0-(SFT_plane(SFT_index+3,4:6))),SFT_plane(SFT_index+3,1:3))<SFTol))
    
      check=1;
      break;
    else
        check=0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

