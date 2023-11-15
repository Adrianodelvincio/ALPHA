{

  // below I look for "243 List 6 Freq 0"
  
  int irun= 68465; 
  double veto_before_time= 0; 

  TA2Plot plotL6(true);  plotL6.AddDumpGates(irun, {"243 List 6 Freq 0"}, {-1});   plotL6.LoadData();
  plotL6.ExportCSV(TString::Format("r%d_laser6_exp_freq0",irun).Data());

  TA2Plot plotMW1(true);  plotMW1.AddDumpGates(irun, {"MW Frq 1 Dump"}, {-1});   plotMW1.LoadData();
  plotMW1.ExportCSV(TString::Format("r%d_uw_exp_freq1",irun).Data());
  
  TA2Plot plotMW2(true);  plotMW2.AddDumpGates(irun, {"MW Frq 2 Dump"}, {-1});   plotMW2.LoadData();
  plotMW2.ExportCSV(TString::Format("r%d_uw_exp_freq2",irun).Data());
  
  TA2Plot plotMW3(true);  plotMW3.AddDumpGates(irun, {"MW Frq 3 Dump"}, {-1});   plotMW3.LoadData();
  plotMW3.ExportCSV(TString::Format("r%d_uw_exp_freq3",irun).Data());
  
  TA2Plot plotMW4(true);  plotMW4.AddDumpGates(irun, {"MW Frq 4 Dump"}, {-1});   plotMW4.LoadData();
  plotMW4.ExportCSV(TString::Format("r%d_uw_exp_freq4",irun).Data());
  
  TA2Plot plotMW5(true);  plotMW5.AddDumpGates(irun, {"MW Frq 5 Dump"}, {-1});   plotMW5.LoadData();
  plotMW5.ExportCSV(TString::Format("r%d_uw_exp_freq5",irun).Data());

  TA2Plot plotMW6(true);  plotMW6.AddDumpGates(irun, {"MW Frq 6 Dump"}, {-1});   plotMW6.LoadData();
  plotMW6.ExportCSV(TString::Format("r%d_uw_exp_freq6",irun).Data());
  
  TA2Plot plotMW7(true);  plotMW7.AddDumpGates(irun, {"MW Frq 7 Dump"}, {-1});   plotMW7.LoadData();
  plotMW7.ExportCSV(TString::Format("r%d_uw_exp_freq7",irun).Data());

  TA2Plot plotMW8(true);  plotMW8.AddDumpGates(irun, {"MW OffResonance Dump"}, {-1});   plotMW8.LoadData();
  plotMW8.ExportCSV(TString::Format("r%d_uw_exp_offres",irun).Data());

  TA2Plot plotMWoff(true);  plotMWoff.AddDumpGates(irun, {"MW Power Off"}, {-1});   plotMWoff.LoadData();
  plotMWoff.ExportCSV(TString::Format("r%d_uw_off",irun).Data());


  TA2Plot plotL5(true);  plotL5.AddDumpGates(irun, {"243 List 5 Freq 0"}, {-1});   plotL5.LoadData();
  plotL5.ExportCSV(TString::Format("r%d_laser5_exp_freq0",irun).Data());

  TA2Plot plotcMW1(true);  plotcMW1.AddDumpGates(irun, {"c MW Frq 1 Dump"}, {-1});   plotcMW1.LoadData();
  plotcMW1.ExportCSV(TString::Format("r%d_c_uw_exp_freq1",irun).Data());
  
  TA2Plot plotcMW2(true);  plotcMW2.AddDumpGates(irun, {"c MW Frq 2 Dump"}, {-1});   plotcMW2.LoadData();
  plotcMW2.ExportCSV(TString::Format("r%d_c_uw_exp_freq2",irun).Data());
  
  TA2Plot plotcMW3(true);  plotcMW3.AddDumpGates(irun, {"c MW Frq 3 Dump"}, {-1});   plotcMW3.LoadData();
  plotcMW3.ExportCSV(TString::Format("r%d_c_uw_exp_freq3",irun).Data());
  
  TA2Plot plotcMW4(true);  plotcMW4.AddDumpGates(irun, {"c MW Frq 4 Dump"}, {-1});   plotcMW4.LoadData();
  plotcMW4.ExportCSV(TString::Format("r%d_c_uw_exp_freq4",irun).Data());
  
  TA2Plot plotcMW5(true);  plotcMW5.AddDumpGates(irun, {"c MW Frq 5 Dump"}, {-1});   plotcMW5.LoadData();
  plotcMW5.ExportCSV(TString::Format("r%d_c_uw_exp_freq5",irun).Data());

  TA2Plot plotcMW6(true);  plotcMW6.AddDumpGates(irun, {"c MW Frq 6 Dump"}, {-1});   plotcMW6.LoadData();
  plotcMW6.ExportCSV(TString::Format("r%d_c_uw_exp_freq6",irun).Data());
  
  TA2Plot plotcMW7(true);  plotcMW7.AddDumpGates(irun, {"c MW Frq 7 Dump"}, {-1});   plotcMW7.LoadData();
  plotcMW7.ExportCSV(TString::Format("r%d_c_uw_exp_freq7",irun).Data());

  TA2Plot plotcMW8(true);  plotcMW8.AddDumpGates(irun, {"c MW OffResonance Dump"}, {-1});   plotcMW8.LoadData();
  plotcMW8.ExportCSV(TString::Format("r%d_c_uw_exp_offres",irun).Data());

  TA2Plot plotcMWoff(true);  plotcMWoff.AddDumpGates(irun, {"c MW Power Off"}, {-1});   plotcMWoff.LoadData();
  plotcMWoff.ExportCSV(TString::Format("r%d_c_uw_off",irun).Data());

  
  TA2Plot plotUWpre0(true);  plotUWpre0.AddDumpGates(irun, {"Microwave Dump"}, {0});   plotUWpre0.LoadData();
  plotUWpre0.ExportCSV(TString::Format("r%d_prelaser_uw_0",irun).Data());
  
  TA2Plot plotUW0(true);  plotUW0.AddDumpGates(irun, {"Microwave Dump"}, {1});   plotUW0.LoadData();
  plotUW0.ExportCSV(TString::Format("r%d_postlaser_uw_0",irun).Data());
  
  TA2Plot plotUW1(true);  plotUW1.AddDumpGates(irun, {"Microwave Dump"}, {2});   plotUW1.LoadData();
  plotUW1.ExportCSV(TString::Format("r%d_postlaser_uw_1",irun).Data()); 



  
  
  // get all shutter-open / shutter-closed pairs
  std::vector<TA2Spill> spillsMix = Get_A2_Spills(irun, {"Mixing"}, {-1});    
  std::vector<TA2Spill> spillsFRD = Get_A2_Spills(irun, {"FastRampDown"}, {0});
  std::vector<TA2Spill> spillsUW = Get_A2_Spills(irun, {"Microwave Dump"}, {-1});    
  
  double tstart = spillsMix.back().GetStopTime();
  double tstop = spillsFRD.back().GetStartTime();
  
  auto shopen = GetSISTimeAndCounts(irun,"LASER_SHUTTER_OPEN",{tstart},{tstop}); 
  auto shclosed = GetSISTimeAndCounts(irun,"LASER_SHUTTER_CLOSE",{tstart},{tstop});   


  double tcyclestart = 604800; // 1 week
  double tcyclestop = -1;

  TA2Plot plot8(true);
  TA2Plot plot9(true);

    
  //  for (int i=0; i<shopen.size(); i++){
  for (int i=shopen.size()-1 ; i>=0; i--){
    
    // now look for spills within a shutter-open / shutter-closed pair ; find the first one and the last one with a laser
    auto s = Get_A2_Spills(irun,shopen.at(i).first,shclosed.at(i).first); 


    bool has_laser = false;
    for (auto k: s){
      if (k.GetSanitisedName() == "243_List_6_Freq_0"){
	has_laser = true;
	break;
      }
      else if (k.GetSanitisedName() == "243_List_5_Freq_0"){
	has_laser = true;
	break;
      }
    }
    
    // inter shutter 
    if (has_laser) {

      if (tcyclestart > shopen.at(i).first)
	tcyclestart = shopen.at(i).first;
      
      if (tcyclestop < shclosed.at(i).first)
	tcyclestop = shclosed.at(i).first;

    }
    else {
      shopen.erase (shopen.begin() +i);
      shclosed.erase (shclosed.begin() +i);
    }

  }

  if (tcyclestart<veto_before_time)
    tcyclestart=veto_before_time;

  if (tcyclestop<veto_before_time)
    tcyclestop=veto_before_time;


  TA2Plot plot11(true);   plot11.AddTimeGate(irun, spillsUW.at(0).GetStopTime(),tcyclestart);   plot11.LoadData();
  plot11.ExportCSV(TString::Format("r%d_postuw_prelaser",irun).Data());                                                                                                                                     
  
  for (int i=0; i<shopen.size(); i++){
    
    if (shopen.at(i).first > tcyclestop  || shclosed.at(i).first < tcyclestart)     // these are all within the same cycle
      continue;

    //    cout << tcyclestart << " " << shopen.at(i).first << " - " << shclosed.at(i).first << " " << tcyclestop << endl; 
    
    if (i>0)
      plot8.AddTimeGate(irun, shclosed.at(i-1).first, shopen.at(i).first); 
    
    auto s = Get_A2_Spills(irun,shopen.at(i).first,shclosed.at(i).first); 

    if (s.size()>0){
   
      std::vector<double> starttimes;
      std::vector<double> stoptimes;
      
      for (auto k: s){
	starttimes.push_back(k.GetStartTime()); 
	stoptimes.push_back(k.GetStopTime()); 
      }  
      std::sort(starttimes.begin(),starttimes.end());
      std::sort(stoptimes.begin(),stoptimes.end());
            
      plot9.AddTimeGate(irun, shopen.at(i).first, starttimes.front());  // shutter open to first startdump
      plot9.AddTimeGate(irun, stoptimes.back(), shclosed.at(i).first ); // last stopdump to shutter closed
      
      for (int j=0; j<starttimes.size()-1; j++) 
        plot9.AddTimeGate(irun, stoptimes.at(j), starttimes.at(j+1) );   // between first dump stop and last dump start
      
    }

  }
  plot8.LoadData();
  plot8.ExportCSV(TString::Format("r%d_shutter_closed",irun).Data()); 

  plot9.LoadData();
  plot9.ExportCSV(TString::Format("r%d_interdumps",irun).Data()); 


  TA2Plot plot6(true);   plot6.AddTimeGate(irun, tcyclestop, spillsUW.at(1).GetStartTime());   plot6.LoadData();
  plot6.ExportCSV(TString::Format("r%d_postlaser_preuw",irun).Data()); 

  TA2Plot plot7(true);   plot7.AddTimeGate(irun, spillsUW.at(1).GetStopTime(), spillsUW.at(2).GetStartTime());   plot7.LoadData();
  plot7.ExportCSV(TString::Format("r%d_inter_uw",irun).Data()); 

  TA2Plot plot5(true);   plot5.AddTimeGate(irun, spillsUW.at(2).GetStopTime(), spillsFRD.at(0).GetStartTime());   plot5.LoadData();
  plot5.ExportCSV(TString::Format("r%d_postuw_prefrd",irun).Data()); 

  
  TA2Plot plot10(true);  plot10.AddDumpGates(irun, {"FastRampDown"}, {0});   plot10.LoadData();
  plot10.ExportCSV(TString::Format("r%d_frd",irun).Data());

   

} 
