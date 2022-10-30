/* Read inputs in input.csv */

#include "fstr.h"

bool Fstr::cookInput(std::string path)
{
  // Creating an object of CSVReader
  CSVReader reader(path);
  bool successful;
  // Get the data from CSV File
  std::vector<std::vector<std::string> > dataList;
  successful = reader.getData(dataList); // All data are stored as string firstly
  if (! successful )
  {
    //std::cout<<"There are something wrong while reading file."<<std::endl;
    //std::cout<<"I will use default parameters to caculate."<<std::endl;
    return false;
  }
  // Split all of the data
  for(auto itRow = dataList.begin(); itRow != dataList.end(); ++itRow)
  {
    auto size=itRow->size();
    //std::cout<<size<<std::endl;
    if(size==1)
    {
      continue;
    }
    auto itCol = itRow->begin() + 1;
    if((*itRow)[0]=="maxL")
    {
      maxL.clear();
      if(std::isdigit((*itCol)[0]) || ((*itCol)[0] == '-' && std::isdigit((*itCol)[1])))
      {
        maxL.push_back(std::stod(*itCol));
        maxL_j=true;
      }
      else
      {
        maxL_j=false;
      }
    }
    if((*itRow)[0]=="maxOpT")
    {
      maxOpT.clear();
      if(std::isdigit((*itCol)[0]) || ((*itCol)[0] == '-' && std::isdigit((*itCol)[1])))
      {
        maxOpT.push_back(std::stod(*itCol));
        maxOpT_j=true;;
      }
      else
      {
        maxOpT_j=false;
      }


    }
    if((*itRow)[0]=="opJ")
    {
      opJ.clear();
      if(std::isdigit((*itCol)[0]) || ((*itCol)[0] == '-' && std::isdigit((*itCol)[1])))
      {
        opJ.push_back(std::stod(*itCol));
        opJ_j=true;
      }
      else
      {
        opJ_j=false;
      }
    }
    if((*itRow)[0]=="parity")
    {
      parity.clear();
      if(std::isdigit((*itCol)[0]) || ((*itCol)[0] == '-' && std::isdigit((*itCol)[1])))
      {
        parity.push_back(std::stod(*itCol));
        parity_j=true;
      }
      else
      {
        parity_j=false;
      }
    }
    if((*itRow)[0]=="label")
    {
      channels.clear();
      for(auto itCol = itRow->begin() + 1; itCol != itRow->end(); ++itCol)
      {
        if(std::isdigit((*itCol)[0]) || ((*itCol)[0] == '-' && std::isdigit((*itCol)[1])))
        {
          channels.push_back(std::stod(*itCol));
          channels_j=true;
        }
        else
        {
          channels_j=false;
        }
      }
    }
    if((*itRow)[0]=="mambda")
    {
      mambda.clear();
      if(std::isdigit((*itCol)[0]) || ((*itCol)[0] == '-' && std::isdigit((*itCol)[1])))
      {
        mambda.push_back(std::stod(*itCol));
        mambda_j=true;
      }
      else
      {
        mambda_j=false;
      }
    }
    if((*itRow)[0]=="NP")
    {
      NP.clear();
      for(auto itCol = itRow->begin() + 1; itCol != itRow->end(); ++itCol)
      {
        if(std::isdigit((*itCol)[0]) || ((*itCol)[0] == '-' && std::isdigit((*itCol)[1])))
        {
          NP.push_back(std::stod(*itCol));
          NP_j=true;
        }
        else
        {
          NP_j=false;
        }
      }
    }
    if((*itRow)[0]=="potential")
    {
      potential.clear();
      potential.push_back(*itCol);
      potential_j=true;

    }
    if((*itRow)[0]=="NX")
    {
      NX.clear();
      if(std::isdigit((*itCol)[0]) || ((*itCol)[0] == '-' && std::isdigit((*itCol)[1])))
      {
        NX.push_back(std::stod(*itCol));
        NX_j=true;
      }
      else
      {
        NX_j=false;
      }
    }
    if((*itRow)[0]=="CMRatio")
    {
      CMRatio.clear();
      if(std::isdigit((*itCol)[0]) || ((*itCol)[0] == '-' && std::isdigit((*itCol)[1])))
      {
        CMRatio.push_back(std::stod(*itCol));
        CMRatio_j=true;
      }
      else
      {
        CMRatio_j=false;
      }
    }
    if((*itRow)[0]=="RFMTHD")
    {
      RFMTHD.clear();
      RFMTHD.push_back(*itCol);
      RFMTHD_j=true;
    }
    if((*itRow)[0]=="BE01")
    {
      BE01.clear();
      for(auto itCol = itRow->begin() + 1; itCol != itRow->end(); ++itCol) {
        if(std::isdigit((*itCol)[0]) || ((*itCol)[0] == '-' && std::isdigit((*itCol)[1])))
        {
          BE01.push_back(std::stod(*itCol));
          BE01_j=true;
        }
        else
        {
          BE01_j=false;
        }
      }
    }
    if((*itRow)[0]=="ErrorBE")
    {
      error.clear();
      if(std::isdigit((*itCol)[0]) || ((*itCol)[0] == '-' && std::isdigit((*itCol)[1])))
      {
        error.push_back(std::stod(*itCol));
        error_j=true;
      }
      else
      {
        error_j=false;
      }
    }
    if((*itRow)[0]=="phi")
    {
      phi.clear();
      if(std::isdigit((*itCol)[0]) || ((*itCol)[0] == '-' && std::isdigit((*itCol)[1])))
      {
        phi.push_back(std::stod(*itCol));
        phi_j=true;
      }
      else
      {
        phi_j=false;
      }
    }
    if((*itRow)[0]=="ratioNQP")
    {
      ratioNQP.clear();
      if(std::isdigit((*itCol)[0]) || ((*itCol)[0] == '-' && std::isdigit((*itCol)[1])))
      {
        ratioNQP.push_back(std::stod(*itCol));
        ratioNQP_j=true;
      }
      else
      {
        ratioNQP_j=false;
      }
    }

  }
  return true;
};

void Fstr::PrintParameter()
{
  if(maxL_j)
  {
    std::cout << "maxL = " << maxL[0] << std::endl;
  }
  else
  {
    std::cout << "maxL = " << maxL[0] <<"(default)"<< std::endl;
  }

  if(maxOpT_j)
  {
    std::cout << "maxOpT = " << maxOpT[0] << std::endl;
  }
  else
  {
    std::cout << "maxOpT = " << maxOpT[0] <<"(default)"<< std::endl;
  }

  if(opJ_j)
  {
    std::cout << "opJ = " << opJ[0] << std::endl;
  }
  else
  {
    std::cout << "opJ = " << opJ[0] <<"(default)"<< std::endl;
  }

  if(parity_j)
  {
    std::cout << "parity = " << parity[0] << std::endl;
  }
  else
  {
    std::cout << "parity = " << parity[0] <<"(default)"<< std::endl;
  }

  if(potential_j)
  {
    std::cout << "Potential = " << potential[0] << std::endl;
  }
  else
  {
    std::cout << "Potential = " << potential[0] <<"(default)"<< std::endl;
  }

  if(mambda_j)
  {
    std::cout << "Mambda = " << mambda[0] << std::endl;
  }
  else
  {
    std::cout << "Mambda = " << mambda[0] <<"(default)"<< std::endl;
  }

  if(NP_j)
  {
    std::cout << "dNP = " << int(NP[0]) << ", NP0 = " << int(NP[1]) << ", NP1 = " << int(NP[2]) << std::endl;
  }
  else
  {
    std::cout << "dNP = " << int(NP[0]) << ", NP0 = " << int(NP[1]) << ", NP1 = " << int(NP[2]) <<"(default)"<< std::endl;
  }

  if(NX_j)
  {
    std::cout << "NX = " << NX[0] << std::endl;
  }
  else
  {
    std::cout << "NX = " << NX[0] <<"(default)"<< std::endl;
  }

  if(CMRatio_j)
  {
    std::cout << "CTanMesh/Mambda = " << CMRatio[0] << std::endl;
  }
  else
  {
    std::cout << "CTanMesh/Mambda = " << CMRatio[0] <<"(default)"<< std::endl;
  }

  if(RFMTHD_j)
  {
    if(RFMTHD[0]=="Scant"||RFMTHD[0]=="Newton")
    {
      std::cout << "RootFind Method = " << RFMTHD[0] << std::endl;
    }
    else
    {
      std::cout <<RFMTHD[0] <<" not recognized, I will use default RootFind Method"<<std::endl;
      RFMTHD[0]="Scant";
    }
  }
  else
  {
    std::cout << "RootFind Method = " << RFMTHD[0] <<"(default)"<< std::endl;
  }

  if(BE01_j)
  {
    std::cout << "Trial values for BE = " << BE01[0]
    << ",  " << BE01[1] << std::endl;
  }
  else
  {
    std::cout << "Trial values for BE = " << BE01[0]
    << ",  " << BE01[1] <<"(default)"<< std::endl;
  }

  if(error_j)
  {
    std::cout << "Error of BE = " << error[0] << std::endl;
  }
  else
  {
    std::cout << "Error of BE = " << error[0] <<"(default)"<< std::endl;
  }

  if ( potential[0].find("cmplx") != std::string::npos )
  {
    if(phi_j)
    {
      std::cout << "phi (degrees) = " << phi[0] << std::endl;
    }
    else
    {
      std::cout << "phi (degrees) = " << phi[0] <<"(default)"<< std::endl;
    }
  }

  if(ratioNQP_j)
  {
    std::cout << "NQ/NP = " << ratioNQP[0] << std::endl;
  }
  else
  {
    std::cout << "NQ/NP = " << ratioNQP[0] <<"(default)"<< std::endl;
  }
}
