// ROOT macrot to make dummy ped and cm files for testing code.

int  make_dummy_ped_and_cm(const std::string pedfilename = "ft_ped.dat", const std::string cmfilename = "ft_cm.dat"){

  const int NPlanes_ft {28};
  std::vector<int> nstrips_ft_byplane = {
    3968, 3456, 3968, 3456, 3840, 3840, 3840, 3840, 3840, 3840, 3840, 3840, 1280, 1536, 1280, 1536, 1280, 1536, 1280, 1536, 1280, 1536, 1280, 1536, 1280, 1536, 1280, 1536
  };

  if (NPlanes_ft!=nstrips_ft_byplane.size()) return(-1);

  std::ofstream pedfile(pedfilename);
  std::ofstream cmfile(cmfilename);
  
  for (int iplane=0; iplane<NPlanes_ft; iplane++){
    pedfile << "GEM_plane# " << iplane << std::endl;
    cmfile << "GEM_plane# " << iplane << std::endl;
    
    for (int istrip=0; istrip<nstrips_ft_byplane.at(iplane); istrip++){
      pedfile << istrip << " " << "20" << " " << "10" << std::endl;
      
      if (istrip%128 == 0){
	int iAPV = istrip/128;
	cmfile << iAPV << " " << "1000" << " " << "20" << std::endl;
      }
    }

    pedfile << std::endl;
    cmfile << std::endl;
  }

  pedfile.close();

  return 0;
}
