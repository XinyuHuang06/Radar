function [rho_0, rho_1] = Update_rho(DataSetPackets, ParameterPackets)
    rho_0 = DataSetPackets.packets.rho_0*2;
    rho_1 = DataSetPackets.packets.rho_1*2;
end