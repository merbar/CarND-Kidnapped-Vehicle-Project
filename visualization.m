pfparticles = importdata("pf_particles.txt");
gtdata = importdata("data/gt_data.txt");

ar_size = size(pfparticles);
frames = ar_size(1) / 20;
pfparticles_perFrame = {};
for i = 1:frames
    particles_per_frame = [];
    for j = 1:20
        index = (i-1)*20+j;
        particle_x = pfparticles(index,1);
        particle_y = pfparticles(index,2);
        particle_theta = pfparticles(index,3);
        particles_per_frame=[particles_per_frame; [particle_x, particle_y, particle_theta]];
    end
    pfparticles_perFrame = [pfparticles_perFrame, particles_per_frame];
end

world_xmin = -70;
world_xmax = 350;
world_ymin = -120;
world_ymax = 60;

clear particles_per_frame index i j particle_x particle_y particle_theta ar_size pfparticles;

%bla = cell2mat(pfparticles_perFrame(1)