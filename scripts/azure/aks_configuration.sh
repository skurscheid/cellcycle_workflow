#!/bin/sh

export myRG="CellCycle"
export myAKS="aks-cellcycle"

{
  "appId": "b14db6fe-8f7a-48f7-a9a2-e8516bb474a6",
  "displayName": "azure-cli-2019-02-13-06-28-46",
  "name": "http://azure-cli-2019-02-13-06-28-46",
  "password": "64dff45e-0513-4924-aa27-c2d4055b3895",
  "tenant": "e37d725c-ab5c-4624-9ae5-f0533e486437"
}

# $az acr show --resource-group "CellCycle" --name acrcellcycle --query "id" --output tsv

/subscriptions/012423da-b029-47f5-9a5a-46caecca2eb2/resourceGroups/CellCycle/providers/Microsoft.ContainerRegistry/registries/acrCellCycle

# $az role assignment create --assignee b14db6fe-8f7a-48f7-a9a2-e8516bb474a6 --scope /subscriptions/012423da-b029-47f5-9a5a-46caecca2eb2/resourceGroups/CellCycle/providers/Microsoft.ContainerRegistry/registries/acrCellCycle --role acrpull

{
  "canDelegate": null,
  "id": "/subscriptions/012423da-b029-47f5-9a5a-46caecca2eb2/resourceGroups/CellCycle/providers/Microsoft.ContainerRegistry/registries/acrCellCycle/providers/Microsoft.Authorization/roleAssignments/acf87719-c163-4d70-b9dd-1ef57324f6fe",
  "name": "acf87719-c163-4d70-b9dd-1ef57324f6fe",
  "principalId": "2f1e4223-3cf6-4e4d-846e-00d5e4520caa",
  "resourceGroup": "CellCycle",
  "roleDefinitionId": "/subscriptions/012423da-b029-47f5-9a5a-46caecca2eb2/providers/Microsoft.Authorization/roleDefinitions/7f951dda-4ed3-4680-a7ca-43fe172d538d",
  "scope": "/subscriptions/012423da-b029-47f5-9a5a-46caecca2eb2/resourceGroups/CellCycle/providers/Microsoft.ContainerRegistry/registries/acrCellCycle",
  "type": "Microsoft.Authorization/roleAssignments"
}

az aks create \
    --resource-group "CellCycle" \
    --name aks-cellcycle \
    --node-count 3 \
    --service-principal b14db6fe-8f7a-48f7-a9a2-e8516bb474a6 \
    --client-secret 64dff45e-0513-4924-aa27-c2d4055b3895 \
    --generate-ssh-keys

az aks get-credentials --resource-group $myRG --name $myAKS

kubectl get nodes

az acr list --resource-group $myRG --query "[].{acrLoginServer:loginServer}" --output table

# create static file share for aks
export AKS_PERS_STORAGE_ACCOUNT_NAME="cellcycledata"
export STORAGE_ACCOUNT_KEY=$(az storage account keys list -n $AKS_PERS_STORAGE_ACCOUNT_NAME -g $myRG --query='[0].value' | tr -d '"')
kubectl create secret generic azure-secret --from-literal=azurestorageaccountname=$AKS_PERS_STORAGE_ACCOUNT_NAME --from-literal=azurestorageaccountkey=$STORAGE_ACCOUNT_KEY
kubectl apply -f azure-files-pod.yaml
kubectl describe pod mypod