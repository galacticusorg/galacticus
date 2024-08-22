# Provides functions used for validation of models.
# Andrew Benson (12-August-2024)
import os
import numpy as np
import h5py
import re
import json
import codecs
from git import Repo

def extract(fileName, name, suffix, parameterFileName):
    # Extract a likelihood measure from a given model for all analyses in that model. The likelihood measure we choose is -logℒ,
    # where we assume that ℒ is "unnormalized" - that is, for a normal distribution it is
    #
    #  ∏ᵢ exp(-½(x-xᵢ)²/σᵢ²)
    #
    # and not
    #
    #  ∑ᵢ exp(-½(x-xᵢ)²/σᵢ²)/√[2πσᵢ²]
    #
    # Then,
    #
    #  -logℒ = ∑ᵢ ½(x-xᵢ)²/σᵢ²
    #
    # and the derivative of this with respect to model predictions is
    #
    #  ∂(-logℒ)/∂x = (x-xᵢ)/σᵢ².
    #
    # Then a fractional change in the offset of the model from the data of Δ(x-xᵢ)/(x-xᵢ) results in a change in our metric of
    #
    #  Δ(-logℒ) = ∑ᵢ (x-xᵢ)²/σᵢ² Δ(x-xᵢ)/(x-xᵢ)
    #
    # If that fractional change is the same, Δf, for all points, i, then
    #
    #  Δ(-logℒ) = Δf ∑ᵢ (x-xᵢ)²/σᵢ² = 2 Δf (-logℒ),
    #
    # or
    #
    #  Δ(-logℒ)/(-logℒ) = Δf ∑ᵢ (x-xᵢ)²/σᵢ² = 2 Δf.
    #
    # Therefore a fractional shift in the model relative to the data results in a corresponding fractional shift in our
    # metric. This allows us to use a percentage change threshold in this metric as a warning for a significant shift in the
    # model even when the model is not a good match to the data (i.e. when the likelihood is low and can change hugely due to
    # even small shifts in the model).
    likelihoods = []
    results     = []
    model       = h5py.File(fileName,"r")
    analyses    = model['analyses']
    for analysisName, analysisGroup in analyses.items():
        attributes = {
            "name"       : analysisName,
            "xAxisLabel" : "x",
            "yAxisLabel" : "y",
            "xAxisIsLog" :  0 ,
            "yAxisIsLog" :  1 ,
            "targetLabel": "data"
        }
        for attributeName in analysisGroup.attrs:
            attributes[attributeName] = analysisGroup.attrs[attributeName]
            if isinstance(attributes[attributeName], bytes):
                attributes[attributeName] = attributes[attributeName].decode('utf-8')
        logLikelihood = attributes['logLikelihood']
        print(analysisName+"\t"+str(logLikelihood))
        likelihoods.append(
            {
        	"name" : name+" - Likelihood - "+analysisName,
        	"unit" : "-logℒ"                            ,
        	"value": str(np.abs(logLikelihood))
            }
            )
        # Skip cases for which we have no "type" specified.
        if not "type" in attributes:
            print("Warning: analysis '"+analysisName+"' has no 'type' attribute, so it can not be processed.")
        elif attributes['type'] == "function1D":
            # Simple 1D function - will be shown as an x-y scatter plot.
            # Validate attributes.
            datasetNames = (
        	{"name": 'xDataset'         , "required": 1},
        	{"name": 'yDataset'         , "required": 1},
        	{"name": 'yDatasetTarget'   , "required": 0},
        	{"name": 'yCovariance'      , "required": 0},
        	{"name": 'yCovarianceTarget', "required": 0},
        	{"name": 'yErrorLower'      , "required": 0},
        	{"name": 'yErrorUpper'      , "required": 0},
        	{"name": 'yErrorLowerTarget', "required": 0},
        	{"name": 'yErrorUpperTarget', "required": 0},
            )
            for dataset in datasetNames:
                if dataset['required'] and not dataset['name'] in attributes:
                    print("Error: attribute '"+dataset['name']+"' is missing from analysis '"+analysisName+"' but is required.")
                    os.abort()
            # Read the datasets.
            data = {}
            for dataset in datasetNames:
                if dataset['name'] in attributes:
                    analysisDatasetName = attributes[dataset['name']]
                    if not analysisDatasetName in analysisGroup.keys():
                        print("Analysis: '"    +analysisName       +"'")
                        print("Generic name: '"+dataset['name']    +"'")
                        print("Actual name: '" +analysisDatasetName+"'")
                        print("Available datasets:")
                        print("\n".join(analysisGroup.keys()))
                        print("failed to find dataset")
                        os.abort()
                    data[dataset['name']] = analysisGroup[analysisDatasetName][:]
            # Extract errors.
            if "yCovariance" in data:
                data['yError'      ] = np.sqrt(np.diagonal(data['yCovariance'      ]))
                del data['yCovariance'      ]
            if "yCovarianceTarget" in data:
                data['yErrorTarget'] = np.sqrt(np.diagonal(data['yCovarianceTarget']))
                del data['yCovarianceTarget']
            # Store results.
            result = {
                "attributes": {},
                "data"      : {}
            }
            for attributeName in attributes:
                if re.search(r"^[xy]AxisIsLog$",attributeName) or attributeName == "logLikelihood":
                    result['attributes'][attributeName] = str(attributes[attributeName])
                else:
                    # LaTeX conversions.
                    attributes[attributeName] = re.sub(r"\$",""                        ,attributes[attributeName])
                    attributes[attributeName] = re.sub(r"\\mathrm\{([^\}]+)\}",r"\g<1>",attributes[attributeName])
                    attributes[attributeName] = re.sub(r"\\hbox\{([^\}]+)\}",r"\g<1>"  ,attributes[attributeName])
                    attributes[attributeName] = re.sub(r"\\odot","☉"                   ,attributes[attributeName])
                    attributes[attributeName] = re.sub(r"\\langle","⟨"                 ,attributes[attributeName])
                    attributes[attributeName] = re.sub(r"\\rangle","⟩"                 ,attributes[attributeName])
                    attributes[attributeName] = re.sub(r"\\star","★"                   ,attributes[attributeName])
                    attributes[attributeName] = re.sub(r"\\log_\{10\}","log₁₀"         ,attributes[attributeName])
                    attributes[attributeName] = re.sub(r"\\sigma","σ"                  ,attributes[attributeName])
                    attributes[attributeName] = re.sub(r"\^\{-1\}","⁻¹"                ,attributes[attributeName])
                    attributes[attributeName] = re.sub(r"\\,"," "                      ,attributes[attributeName])
                    result['attributes'][attributeName] = attributes[attributeName]
            for dataName in data:
                result['data'][dataName] = list(data[dataName])
            results.append(result)
    # Write benchmark results.
    f = codecs.open("outputs/validate_"+suffix+".json", "w", "utf-8")
    f.write(json.dumps(likelihoods,indent=4,ensure_ascii=False))
    f.close()
    # Interface with git.
    repo         = Repo(os.environ['GALACTICUS_EXEC_PATH'])
    actor        = repo.head.commit.author
    lastRevision = repo.head.object.hexsha
    authorName   = actor.name
    authorEmail  = actor.email
    authorDate   = str(repo.head.commit.committed_datetime)
    message      = repo.head.commit.message
    # Write results.
    output = {
        "repoUrl"      : "https://github.com/galacticusorg/galacticus",
        "parameterFile": parameterFileName,
        "commit"       : {
            "author":  {
                "name" : authorName,
                "email": authorEmail
            },
            "id"       : lastRevision,
            "message"  : message,
            "timestamp": authorDate,
            "url"      : "https://github.com/galacticusorg/galacticus/commit/"+lastRevision
        },
        "results": results
    }
    f = codecs.open("outputs/results_"+suffix+".json", "w", "utf-8")
    f.write("window.ANALYSES_DATA = ")
    f.write(json.dumps(output,indent=4,ensure_ascii=False))
    f.close()
