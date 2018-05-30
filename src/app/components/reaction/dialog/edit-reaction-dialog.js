angular
    .module('app')
    .controller('EditReactionDialogController', EditReactionDialogController);
function EditReactionDialogController($scope, $mdDialog, $http, $log) {
    var self = this;
    self.states = loadAll();
    self.querySearch = querySearch;
    self.selectedItemChange = selectedItemChange;
    self.searchTextChange = searchTextChange;
    $scope.description = "";
    $scope.hide = function () {
        $mdDialog.hide();
    };

    $scope.cancel = function () {
        $mdDialog.cancel();
    };

    $scope.create = function () {
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
            swal('Fail to update reaction','','error')
        })

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
            self.name = item.name;
            self.id = item.id;
            self.description = item.description;
            self.smart = item.smart;
            $log.info('Item changed to ' + JSON.stringify(item));
        }
    }
}