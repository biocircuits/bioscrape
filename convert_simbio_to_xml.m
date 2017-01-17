function []=convert_simbio_to_xml(filename, model)

f = fopen(filename, 'w');
m = model;
fprintf(f, '<model>\n\n');

function [clean_name] = F(dirty_name)
    clean_name = dirty_name;
    clean_name(clean_name == '-') = '_';
    clean_name(clean_name == ' ') = '_';
end

% Go through reactions one at a time. Add the reactions and propensity, but
% do not update the output until later.
number_of_reactions = length(m.Reactions);

param_names = {};
param_values = [];
species_used = {};
for i=1:number_of_reactions,
    r = m.Reactions(i);
    number_of_reactants = length(r.Reactants);
    number_of_products = length(r.Products);
    
    % Create reactants together.
    reactant_string = [];
    for j=1:number_of_reactants,
        reactant_string = [reactant_string F(r.Reactants(j).Name)];
        species_used = [species_used F(r.Reactants(j).Name)];
        if j < length(r.Reactants),
            reactant_string = [reactant_string ' + '];
        end
    end
   
    % Create products together.
    product_string = [];
    for j=1:number_of_products,
        product_string = [product_string F(r.Products(j).Name)];
        species_used = [species_used F(r.Products(j).Name)];
        if j < length(r.Products),
            product_string = [product_string ' + '];
        end
    end
    
    % Get the reaction rate type and parameters.
    kl = r.KineticLaw;
    assert( strcmp(kl.KineticLawName,'MassAction') );
    parameters = r.KineticLaw.getparameters();
    
    % Add parameters to the list of total parameters
    for k=1:length(parameters),
        if ~ismember(parameters(k).Name, param_names),
            param_names = [param_names F(parameters(k).Name)];
            param_values = [param_values parameters(k).Value];
        end
    end
    
    % Forward reaction first.
    
    param = parameters(1).Name;
    fprintf(f, '<reaction text="%s--%s" after="--">\n', ...
               reactant_string, product_string);
    
    species_string = '';     
    for reactant_index = 1:number_of_reactants,
        species_string = [species_string F(r.Reactants(reactant_index).Name)];
        if reactant_index < number_of_reactants,
            species_string = [species_string '*'];
        end
    end
    fprintf(f, '    <propensity type="massaction" k="%s" species="%s" />\n', ...
                   param, species_string);
    fprintf(f,'    <delay type="none" />\n</reaction>\n\n');
    
    % The the reverse reaction, but skip if it does not exist!!
    if length(parameters) < 2,
        continue;
    end
    param = parameters(2).Name;
    
    fprintf(f, '<reaction text="%s--%s" after="--">\n', ...
               product_string, reactant_string);
    
    species_string = '';     
    for product_index = 1:number_of_products,
        species_string = [species_string F(r.Products(product_index).Name)];
        if product_index < number_of_products,
            species_string = [species_string '*'];
        end
    end
    fprintf(f, '    <propensity type="massaction" k="%s" species="%s" />\n', ...
                   param, species_string);
    
    fprintf(f,'    <delay type="none" />\n</reaction>\n\n');
end
        
% Take care of parameters as well.
fprintf(f,'\n\n\n');
for i=1:length(param_values),
    fprintf(f, '<parameter name="%s" value="%f" />\n', ...
            param_names{i}, param_values(i) );
end

% Take care of species last.
fprintf(f,'\n\n\n');
number_of_species = length(m.Species);
for i=1:number_of_species,
   if ismember(F(m.Species(i).Name), species_used),    
      fprintf(f,'<species name="%s" value="%f" />\n', F(m.Species(i).Name), m.Species(i).InitialAmount);
   end
end

fprintf(f, '\n\n</model>\n');

end