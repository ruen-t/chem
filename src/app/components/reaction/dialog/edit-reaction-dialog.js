angular
    .module('app')
    .controller('EditReactionDialogController', EditReactionDialogController);
function EditReactionDialogController($scope, $mdDialog, $http, $log) {
    var self = this;
    self.states = loadAll();
    self.querySearch = querySearch;
    self.selectedItemChange = selectedItemChange;
    self.searchTextChange = searchTextChange;
    self.description = "";
    self.showError = false;
    self.id = -1;
    $scope.hide = function () {
        $mdDialog.hide();
    };
    $scope.cancel = function () {
        $mdDialog.cancel();
    };
    $scope.delete = function () {
        if (self.id < 0) {
            self.showError = true;
        } else {
            $mdDialog.cancel();
            var deleted = false;
            swal({
                title: 'Please type SMARTS of the item which you would like to delete',
                input: 'text',
                inputAttributes: {
                    autocapitalize: 'off'
                },
                showCancelButton: true,
                confirmButtonText: 'Confirm',
                showLoaderOnConfirm: true,
                preConfirm: (value) => {
                    console.log(value)
                    
                },
                allowOutsideClick: false
            }).then((result) => {
                if(result.dismiss){
                  return;  
                }
                if (result.value != self.smart) {
                    swal('Fail to delete reaction', 'Please check your confirm smarts input', 'error')
                } else {
                    $http({
                        url: REACTION_API + "/" + self.id,
                        method: 'DELETE',
                        headers: { 'Content-Type': undefined },
                        transformRequest: angular.identity
                    }).then(function (response) {
                        console.log(response);
                        deleted = true;
                        swal('Successfully delete reaction', '', 'success')
                    }, function (error) {
                        swal('Fail to delete reaction', 'Server error', 'error')
                    })
                }
            }, (error) => {
                swal('Fail to delete reaction', 'Please check your confirm smarts input', 'error')
            })
        }
    };

    $scope.update = function () {
        if (id < 0) {
            self.showError = false;
        } else {
            $mdDialog.hide();
            var name = self.name;
            var smart = self.smart;
            var description = self.description;
            var payload = new FormData();
            var reaction = { name, smart, description };
            payload.append("reaction", JSON.stringify(reaction));

            $http({
                url: REACTION_API + "/" + self.id,
                method: 'POST',
                data: payload,
                headers: { 'Content-Type': undefined },
                transformRequest: angular.identity
            }).then(function (response) {
                console.log(response);
                swal('Successfully update reaction', '', 'success')
            }, function (error) {
                swal('Fail to update reaction', '', 'error')
            })
        }


    };
    function loadAll() {
        return new Promise(function (resolve, reject) {
            $http.get(REACTION_API).then(function (response) {
                self.reactionList = response.data.reaction;
                self.reactionList.forEach(function (element) {
                    element.display = element.name + " " + element.smart;
                });
                resolve()
            });
        })
    }
    function querySearch(query) {
        //console.log(self.reactionList)
        if (typeof self.reactionList == 'undefined') {
            var deferred = $q.defer();
            loadAll().then((success) => {
                var results = query ? self.reactionList.filter(createFilterFor(query)) : self.reactionList;
                deferred.resolve(results)
            })
            return deferred.promise;
        } else {
            var results = query ? self.reactionList.filter(createFilterFor(query)) : self.reactionList;
            return results

        }

    }
    function searchTextChange(text) {
        $log.info('Text changed to ' + text);
    }

    function selectedItemChange(item) {
        if (item != null) {
            self.showError = false;
            self.name = item.name;
            self.id = item.id;
            self.description = item.description;
            self.smart = item.smart;
            $log.info('Item changed to ' + JSON.stringify(item));
        } else {
            self.name = "";
            self.id = -1;
            self.description = "";
            self.smart = "";
        }
    }
}